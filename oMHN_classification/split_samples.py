import fastmhn
import numpy as np

from parameters import *

othetas = [
    np.loadtxt(f"./theta_group{group}.dat") for group in range(nr_groups)
]
data = np.genfromtxt(sample_matrix_path, skip_header=1, dtype=int)
data = data[:, 2:]

sample_Ps = np.zeros((nr_groups, data.shape[0]))

thetas = [
    np.zeros((othetas[0].shape[1], othetas[0].shape[1]))
    for _ in range(nr_groups)
]
clusterings = [None for _ in range(nr_groups)]
for i in range(nr_groups):
    # define thetas and get clusterings
    d = othetas[i].shape[1]
    for e1 in range(d):
        for e2 in range(d):
            if e1 == e2:
                thetas[i][e1, e2] = othetas[i][e1, e2]
            else:
                thetas[i][e1, e2] = othetas[i][e1, e2] - othetas[i][-1, e2]

# loop over all groups
for group in range(nr_groups):
    sample_Ps[group] = np.exp(
        np.array(
            fastmhn.approx.__get_approx_gradient_and_score_contributions(
                thetas[group], data, max_cluster_size=15
            )[1]
        )
    )
    print(f"pTheta entries for group {group} calculated")

sample_Ps = sample_Ps / sample_Ps.sum(axis=0, keepdims=True)

np.savetxt(f"sample_Ps.dat", sample_Ps)
