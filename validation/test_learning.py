import os

import fastmhn
import numpy as np

from parameters import *

os.mkdir("learned_thetas")

max_number_of_events = 24
adam_parameters = {"alpha": 1e-1, "beta1": 0.7, "beta2": 0.9, "N_max": 100}

for nr_theta in range(nr_thetas):
    data = np.loadtxt(f"validation_data/data{nr_theta}.dat")
    data_exact = data[np.sum(data, axis=1) <= max_number_of_events]

    reg = 1 / data_exact.shape[0]

    # exact learning process
    no_clustering = lambda theta: [[i for i in range(d)]]
    theta = fastmhn.learn.learn_mhn(
        data_exact,
        reg=reg,
        gradient_and_score_params={
            "clustering_algorithm": no_clustering,
        },
        adam_params=adam_parameters,
    )
    np.savetxt(f"./learned_thetas/theta{nr_theta}_exact.dat", theta)
    print(f"Exact learning process done")

    # approximate learning processes
    for dmax in dmaxs:
        theta = fastmhn.learn.learn_mhn(
            data_exact,
            reg=reg,
            gradient_and_score_params={
                "max_cluster_size": dmax,
            },
            adam_params=adam_parameters,
        )
        np.savetxt(f"./learned_thetas/theta{nr_theta}_dmax{dmax}.dat", theta)
        print(f"Approximate learning process with dmax={dmax} done")
