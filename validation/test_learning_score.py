import fastmhn
import numpy as np

from parameters import *

max_number_of_events = 24

no_clustering = lambda theta: [[i for i in range(d)]]

for nr_theta in range(nr_thetas):
    data = np.loadtxt(f"validation_data/data{nr_theta}.dat")
    data_exact = data[np.sum(data, axis=1) <= max_number_of_events]

    # exact score
    try:
        theta = np.loadtxt(f"./learned_thetas/theta{nr_theta}_exact.dat")
        score = fastmhn.approx.approx_gradient_and_score(
            theta, data_exact, clustering_algorithm=no_clustering
        )[1]
    except:
        score = np.nan

    print(f"{nr_theta} exact {score}")

    # approximate learning processes
    for dmax in dmaxs:
        try:
            theta = np.loadtxt(
                f"./learned_thetas/theta{nr_theta}_dmax{dmax}.dat"
            )
            score = fastmhn.approx.approx_gradient_and_score(
                theta, data_exact, clustering_algorithm=no_clustering
            )[1]
        except:
            score = np.nan

        print(f"{nr_theta} {dmax} {score}")
