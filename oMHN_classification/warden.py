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


time.sleep(2)

log = dict()
if os.path.exists("log.json"):
    log = json.load(open("log.json", "r"))
else:
    # start a new learning process from scratch
    print("Starting a new learning process from scratch")
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

if log["it"] == nr_iterations:
    print("Learning process finished")
    quit()

if log["phase"] == "CV-start":

    print(f"Initial CV runs for iteration {log['it']}")

    subprocess.run(["python", "./aggregate_CV_runs.py"], check=True)
    log["phase"] = "CV"
    dump_json()

if log["phase"] == "CV":

    print(f"CV for iteration {log['it']} in progress")

    # start CV runs that are still missing
    result = subprocess.run(
        ["python", "./aggregate_CV_runs.py"], capture_output=True, text=True
    )
    time.sleep(1)
    with open("./saved_aggregate.dat", "r") as f:
        aggregate = f.read()
    while color.RED in aggregate:
        result = subprocess.run(
            ["python", "./aggregate_CV_runs.py"],
            capture_output=True,
            text=True,
        )
        with open("./saved_aggregate.dat", "r") as f:
            aggregate = f.read()
    log["phase"] = "Final"
    dump_json()

if log["phase"] == "Final":

    print(
        f"CV for iteration {log['it']} finished, running final "
        "model training for expectation step"
    )

    # get the theta matrices
    while not all(
        [os.path.exists(f"theta_group{g}.dat") for g in range(nr_groups)]
    ):
        subprocess.run(["python", "./get_models.py"], check=True)

    subprocess.run(["python", "./get_models.py"], check=True)

    log["phase"] = "Split"
    dump_json()

if log["phase"] == "Split":

    print(
        f"Model training for iteration {log['it']} done, "
        "performing maximization step"
    )

    # distribute the patients to the samples
    subprocess.run(["rm", "-rf", "sample_Ps.dat"], check=True)
    if slurm:
        subprocess.run(["sbatch", "./srun_split_samples"], check=True)
    else:
        subprocess.run(["python", "-u", "split_samples.py"], check=True)
    while not os.path.exists("./sample_Ps.dat"):
        time.sleep(20)
    log["phase"] = "Done"
    dump_json()


if log["phase"] == "Done":

    print(f"Iteration {log['it']} done, moving to next iteration")

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
