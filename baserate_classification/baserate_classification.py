import fastmhn
import numpy as np

np.random.seed(42)

data = np.genfromtxt(
    "../dataset/sample_matrices/sample_matrix_Pan.txt",
    skip_header=1,
    dtype=int,
)
data = data[:, 2:]

perm = np.random.permutation(data.shape[0])
data = data[perm]

nr_samples = data.shape[0]
d = data.shape[1]

print("| **Clusters** | **AIC**  | **BIC**  |")

for nr_groups in range(2, 51):
    weights = np.ones((nr_groups, nr_samples)) / nr_groups
    freqs = np.zeros((nr_groups, d))
    weights = np.ones((nr_groups, nr_samples)) / nr_groups
    weights += np.random.normal(scale=1e-2, size=weights.shape)
    weights /= np.sum(weights, axis=0, keepdims=True)

    # >>> classification iterations
    Ps = np.ones((nr_groups, nr_samples))
    for classification_iteration in range(20):
        # >>> get the models
        thetas = [
            fastmhn.utility.create_indep_model(data, weights=weights[group])
            for group in range(nr_groups)
        ]
        # <<< get the models

        # >>> get the model probabilities
        logPs = np.array(
            [
                fastmhn.approx.__get_approx_gradient_and_score_contributions(
                    thetas[group], data, max_cluster_size=15
                )[1]
                for group in range(nr_groups)
            ]
        )
        # <<< get the model probabilities

        # >>> get weights
        Ps = np.exp(logPs)
        weights = Ps / np.sum(Ps, axis=0, keepdims=True)
        # <<< get weights

    # <<< classification iterations

    # >>> output AIC
    AIC = 2 * nr_groups * d * d - 2 * np.sum(weights * np.log(Ps))
    BIC = np.log(data.shape[0]) * nr_groups * d * d - 2 * np.sum(
        weights * np.log(Ps)
    )

    print(f"| {str(nr_groups) + ' clusters':<12} | {AIC:.1f} | {BIC:.1f} |")
    # <<< output AIC

    # >>> save final weights
    weights = weights[:, np.argsort(perm)]
    np.savetxt(f"sample_Ps_baserate_{nr_groups}groups.dat", weights)
    # <<< save final weights
