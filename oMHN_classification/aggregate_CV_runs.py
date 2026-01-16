import glob
import json
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from parameters import *


def start_run(group, reg):
    global started

    if slurm:
        subprocess.run(
            [
                "sbatch",
                f"--job-name=CV{group}-{reg}",
                "srun_crossvalidation",
                f"{group}",
                reg,
            ],
            check=True,
        )
    else:
        with open(f"slurm-CV-group{group}-reg{reg}.out", "w") as f:
            subprocess.run(
                [
                    "python",
                    "-u",
                    "./learn_approx_omhn_crossvalidated.py",
                    f"{group}",
                    reg,
                ],
                stdout=f,
                text=True,
            )
    started[group].append(reg)


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


sys.stdout = open("saved_aggregate.dat", "w")

started = [[] for _ in range(nr_groups)]
if os.path.exists("started.json"):
    started = json.load(open("started.json", "r"))

slurm_files = glob.glob("slurm-*")

CV_scores = [dict() for _ in range(nr_groups)]
offsets = [0 for _ in range(nr_groups)]

for file in slurm_files:
    try:
        with open(file, "rb") as f:
            f.seek(-2, 2)
            while f.read(1) != b"\n":
                f.seek(-2, 1)
            last_line = f.readline().decode().strip().split()
        if last_line[0] != "Group":
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
        if "1e-3" not in started[group]:
            start_run(group, "1e-3")
        elif "3e-3" not in started[group]:
            start_run(group, "3e-3")
        elif "1e-2" not in started[group]:
            start_run(group, "1e-2")
        print(
            color.RED
            + f"Group {group} not done, waiting for runs to finish!"
            + color.END
            + f" (score offset {offsets[group]})"
        )
        continue

    best_reg = max(CV_scores[group], key=CV_scores[group].get)
    at_border = (best_reg == max(CV_scores[group])) or (
        best_reg == min(CV_scores[group])
    )
    next_reg = ""

    reg_string = f"{best_reg:.0e}"
    base, exp = reg_string.split("e")
    exp = int(exp)
    if best_reg == min(CV_scores[group]):
        if base == "3":
            next_reg = f"1e{exp}"
        else:
            next_reg = f"3e{exp-1}"
    elif best_reg == max(CV_scores[group]):
        if base == "3":
            next_reg = f"1e{exp+1}"
        else:
            next_reg = f"3e{exp}"

    if next_reg == "":
        if len(CV_scores[group]) == len(started[group]):
            print(
                color.GREEN
                + f"Group {group} done, best reg is "
                + color.BOLD
                + f"{best_reg:.0e}"
                + color.END
                + f" (score offset {offsets[group]})"
            )
        else:
            print(
                color.RED
                + f"Group {group} not done, waiting for runs to finish!"
                + color.END
                + f" (score offset {offsets[group]})"
            )
    else:
        print(
            color.RED
            + f"Group {group} not done, running {next_reg} next!"
            + color.END
            + f" (score offset {offsets[group]})"
        )
        if next_reg not in started[group]:
            start_run(group, next_reg)

    scores = [
        (
            color.BOLD + f"{CV_scores[group][k]:.4f}" + color.END
            if k == best_reg
            else f"{CV_scores[group][k]:.4f}"
        )
        for k in sorted(CV_scores[group])
    ]
    print(f"\tscores: {' '.join(scores)}")

with open("started.json", "w") as f:
    json.dump(started, f)
