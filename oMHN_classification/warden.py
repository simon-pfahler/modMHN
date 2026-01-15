import json
import os
import subprocess
import time
from datetime import datetime

import numpy as np

from parameters import *


def dump_json():
    with open("log.json", "w") as f:
        json.dump(log, f)


def end():
    time.sleep(1)
    dump_json()

    subprocess.run(["bash", "run_classification"])

    exit()


class color:
    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


time.sleep(5)

log = dict()
if os.path.exists("log.json"):
    log = json.load(open("log.json", "r"))
else:
    # start a new learning process from scratch
    log["it"] = 0
    subprocess.run(["mkdir", "sample_Ps"], check=True)
    subprocess.run(["mkdir", "aggregates"], check=True)
    with open(sample_matrix_path, "r") as f:
        lines = f.readlines()

    Ps = np.ones((nr_groups, len(lines) - 1)) / nr_groups
    Ps += np.random.normal(scale=1e-2, size=Ps.shape)
    Ps /= np.sum(Ps, axis=0, keepdims=True)

    np.savetxt("./sample_Ps/sample_Ps_it-1.dat", Ps)
    np.savetxt("./sample_Ps.dat", Ps)
    log["phase"] = "CV-start"
    dump_json()
    time.sleep(10)

if log["it"] == nr_iterations:
    quit()

if log["phase"] == "CV-start":
    # start more runs at once at the start of an iteration
    subprocess.run(["python", "./aggregate_CV_runs.py"], check=True)
    time.sleep(1)

    log["phase"] = "CV"
    dump_json()
    time.sleep(10)

if log["phase"] == "CV":
    # start CV runs that are still missing
    result = subprocess.run(
        ["python", "./aggregate_CV_runs.py"], capture_output=True, text=True
    )
    while color.RED in result.stdout:
        time.sleep(30)
        result = subprocess.run(
            ["python", "./aggregate_CV_runs.py"],
            capture_output=True,
            text=True,
        )
    log["phase"] = "Final"
    dump_json()
    time.sleep(10)

if log["phase"] == "Final":
    # get the theta matrices
    while not all(
        [os.path.exists(f"theta_group{g}.dat") for g in range(nr_groups)]
    ):
        subprocess.run(["python", "./get_models.py"], check=True)
        time.sleep(30)

    subprocess.run(["python", "./get_models.py"], check=True)

    log["phase"] = "Split"
    dump_json()
    time.sleep(10)

if log["phase"] == "Split":
    # distribute the patients to the samples
    subprocess.run(["rm", "-rf", "sample_Ps.dat"], check=True)
    if slurm:
        subprocess.run(["sbatch", "./srun_split_samples"], check=True)
    else:
        subprocess.run(["python", "-u", "split_samples.py"], check=True)
    while not os.path.exists("./sample_Ps.dat"):
        time.sleep(30)
    log["phase"] = "Done"
    dump_json()
    time.sleep(10)


if log["phase"] == "Done":
    # go to next iteration
    subprocess.run(["rm", "-f", "started.json"], check=True)
    subprocess.run(["rm", "-f", "models.json"], check=True)
    subprocess.run(["mkdir", f"thetas_it{log['it']}"], check=True)
    subprocess.run(
        f"mv theta_group*.dat thetas_it{log['it']}", shell=True, check=True
    )
    subprocess.run(
        [
            "mv",
            "saved_aggregate.dat",
            f"aggregates/saved_aggregate_it{log['it']}.dat",
        ],
        check=True,
    )
    subprocess.run(
        ["cp", "sample_Ps.dat", f"sample_Ps/sample_Ps_it{log['it']}.dat"],
        check=True,
    )
    subprocess.run("rm -f slurm-*", shell=True, check=True)
    log["phase"] = "CV-start"
    log["it"] += 1
    end()
