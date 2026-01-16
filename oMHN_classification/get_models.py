import glob
import json
import os
import subprocess
import sys

import numpy as np

from parameters import *


def start_run(group, reg):
    if slurm:
        subprocess.run(
            [
                "sbatch",
                f"--job-name=Model{group}",
                "srun_model",
                f"{group}",
                reg,
            ],
            check=True,
        )
    else:
        with open(f"slurm-Model-group{group}-reg{reg}.out", "w") as f:
            subprocess.run(
                ["python", "-u", "./learn_approx_omhn.py", f"{group}", reg],
                stdout=f,
                text=True,
            )


sys.stdout = open("saved_aggregate.dat", "w")


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


models = [0 for _ in range(nr_groups)]
if os.path.exists("models.json"):
    models = json.load(open("models.json", "r"))

slurm_files = glob.glob("slurm-*")

CV_scores = [dict() for _ in range(nr_groups)]
lls = ["missing" for _ in range(nr_groups)]
offsets = [0 for _ in range(nr_groups)]

for file in slurm_files:
    try:
        with open(file, "rb") as f:
            f.seek(-2, 2)
            while f.read(1) != b"\n":
                f.seek(-2, 1)
            last_line = f.readline().decode().strip().split()
        if last_line[0] != "Group":
            if last_line[0] != "Model":
                continue
            group = int(last_line[1])
            ll = float(last_line[-1])
            lls[group] = ll
            continue
        group = int(last_line[1][:-1])
        reg = float(last_line[-4][:-1])
        score = float(last_line[-3])
        offset = float(last_line[-1][:-1])
        CV_scores[group][reg] = score
        offsets[group] = offset
    except:
        pass

for group in range(nr_groups):
    if len(CV_scores[group]) == 0:
        print(color.RED + f"Group {group} not done!" + color.END)
        continue

    best_reg = max(CV_scores[group], key=CV_scores[group].get)
    at_border = (best_reg == max(CV_scores[group])) or (
        best_reg == min(CV_scores[group])
    )
    next_reg = ""

    if not at_border:
        print(
            color.GREEN
            + f"Group {group} done, best reg is "
            + color.BOLD
            + f"{best_reg:.0e}"
            + color.END
            + f" (score offset {offsets[group]}, log-likelihood {lls[group]})"
        )
        if models[group] == 0:
            start_run(group, f"{best_reg:.0e}")
            models[group] = 1
    else:
        print(
            color.RED
            + f"Group {group} not done!"
            + color.END
            + f" (score offset {offsets[group]}, log-likelihood {lls[group]})"
        )

    scores = [
        (
            color.BOLD + f"{CV_scores[group][k]:.4f}" + color.END
            if k == best_reg
            else f"{CV_scores[group][k]:.4f}"
        )
        for k in sorted(CV_scores[group])
    ]
    print(f"\tscores: {' '.join(scores)}")

try:
    LL = sum(lls)
    params = 47 * 48 * nr_groups
    eff_params = 0
    for nr_group in range(nr_groups):
        theta = np.loadtxt(f"theta_group{nr_group}.dat")
        eff_params += np.sum(np.abs(theta) > 0.05)
    nr_samples = 15816
    with open("ICs.txt", "w") as f:
        f.write(f"AIC = {2*params-2*LL}, effective {2*eff_params-2*LL}\n")
        f.write(
            f"BIC = {params*np.log(nr_samples)-2*LL}, effective {eff_params*np.log(nr_samples)-2*LL}\n"
        )
except:
    print("Creating ICs.txt failed")

with open("models.json", "w") as f:
    json.dump(models, f)
