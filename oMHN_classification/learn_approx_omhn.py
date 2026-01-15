import sys

import fastmhn
import numpy as np

from parameters import *

print(
    f"Starting model learning run for group {sys.argv[1]} and reg {sys.argv[2]}"
)

results_filename = f"theta_group{sys.argv[1]}.dat"
weights_filename = f"sample_Ps.dat"

data = np.genfromtxt(sample_matrix_path, skip_header=1)
data = data[:, 2:]

weights = np.loadtxt(weights_filename)[int(sys.argv[1])]

d = data.shape[1]
N = data.shape[0]

# >>> optimization parameters
gradient_and_score_params = {"max_cluster_size": dmax}
reg = float(sys.argv[2])
adam_params = {
    "alpha": 0.1,
    "beta1": 0.7,
    "beta2": 0.9,
    "eps": 1e-8,
    "N_max": 100,
    "verbose": True,
}
# <<< optimization parameters

# >>> Print information about dataset
avg_MB = np.mean(np.sum(data, axis=1))
max_MB = np.max(np.sum(data, axis=1))
nr_samples_approx = np.sum(
    np.sum(data, axis=1) > gradient_and_score_params["max_cluster_size"]
)
print(
    f"Dataset information:\n"
    f"\t{data.shape[0]} Patients\n"
    f"\tAverage mutational burden: {avg_MB}\n"
    f"\tMaximum mutational burden: {max_MB}\n"
    f"\tNumber of samples with MB > {gradient_and_score_params['max_cluster_size']}: {nr_samples_approx}"
)
# <<< Print information about dataset

theta = fastmhn.learn.learn_omhn(
    data,
    weights=weights,
    reg=reg,
    gradient_and_score_params=gradient_and_score_params,
    adam_params=adam_params,
)

_, score = fastmhn.approx.approx_gradient_and_score(
    fastmhn.utility.cmhn_from_omhn(theta),
    data,
    weights=weights,
    max_cluster_size=dmax,
)

print(f"Model {sys.argv[1]} done, log likelihood is {np.sum(weights)*score}")

np.savetxt(results_filename, theta)
