import os

import fastmhn
import numpy as np

from parameters import *

os.mkdir("logpTheta_entries")
os.mkdir("gradient_entries")

for nr_theta in range(nr_thetas):
    theta_GT = np.loadtxt(f"validation_data/theta{nr_theta}.dat")
    data = np.loadtxt(f"validation_data/data{nr_theta}.dat")
    data = np.unique(data, axis=0)

    # exact logpTheta entries
    gradient_entries_exact = np.zeros((data.shape[0], d, d))
    if not os.path.isfile(f"logpTheta_entries/theta{nr_theta}_exact.npy"):
        logpTheta_entries = np.zeros(data.shape[0])
        logpTheta_entries[:] = np.inf
        gradient_entries_exact[:] = np.inf
        indices_exacts = [
            np.where(data.sum(axis=1) == mb)[0] for mb in list(range(25))
        ]
        for mb in list(range(25)):
            if len(indices_exacts[mb]) == 0:
                continue
            indices_exact = indices_exacts[mb]
            print(
                f"Theta {nr_theta}, exact (mb={mb}): {len(indices_exact)} samples"
            )
            data_exact = data[indices_exact]
            gradients, scores = (
                fastmhn.approx.__get_approx_gradient_and_score_contributions(
                    theta_GT, data_exact, max_cluster_size=None
                )
            )
            logpTheta_entries[indices_exact] = scores
            gradient_entries_exact[indices_exact] = gradients
            print(f"Theta {nr_theta}, mb={mb} done")
            np.save(
                f"logpTheta_entries/theta{nr_theta}_exact", logpTheta_entries
            )

    # approx logpTheta entries
    for dmax in dmaxs:
        if not os.path.isfile(
            f"gradient_entries/theta{nr_theta}_dmax{dmax}.npy"
        ):
            logpTheta_entries = np.zeros(data.shape[0])
            gradient_entries = np.zeros((data.shape[0]))
            print(f"Theta {nr_theta}, dmax {dmax}: {data.shape[0]} samples")
            gradients, scores = (
                fastmhn.approx.__get_approx_gradient_and_score_contributions(
                    theta_GT, data, max_cluster_size=dmax
                )
            )
            logpTheta_entries = scores
            gradient_entries = np.linalg.norm(
                gradient_entries_exact - gradients, axis=(1, 2)
            ) / np.linalg.norm(gradient_entries_exact, axis=(1, 2))
            print(f"Theta {nr_theta}, dmax {dmax} done")
            np.save(
                f"logpTheta_entries/theta{nr_theta}_dmax{dmax}",
                logpTheta_entries,
            )
            np.save(
                f"gradient_entries/theta{nr_theta}_dmax{dmax}", gradient_entries
            )
