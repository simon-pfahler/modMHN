import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# >>> get sample matrix
data = np.genfromtxt(
    "../dataset/sample_matrices/sample_matrix_Pan.txt",
    skip_header=1,
    dtype=int,
)
data = data[:, 2:]
d = data.shape[1]
# <<< get sample matrix

# >>> get gene panel
with open("../dataset/gene_panel.txt") as f:
    eventnames = f.read().strip().split(",")
# <<< get gene panel

groups = 13

group_names_MHN = [f"Group-{i+1}" for i in range(groups)]
group_names_CBN = [f"Group-{i+1}" for i in range(groups)]
group_names_baserate = [f"Group-{i+1}" for i in range(groups)]

# >>> get membership probabilities
Ps_baserate = np.loadtxt(
    f"../classification_results/sample_Ps_baserate_{groups}groups.dat"
)
Ps_MHN = np.loadtxt(
    f"../classification_results/sample_Ps_oMHN_{groups}groups.dat"
)
Ps_CBN = np.loadtxt(
    f"../classification_results/sample_Ps_CBN_{groups}groups.dat"
)
# <<< get membership probabilities

delta_MHN = list()
delta_CBN = list()
delta_baserate = list()


def comp(p1, p2):
    return np.sum(p1 * np.log(p1 / p2))  # KL divergence


for sample_index in tqdm(range(Ps_MHN.shape[1])):
    for e in range(d):
        test_sample = data[sample_index]
        if test_sample[e] == 0:
            continue
        predecessor = test_sample.copy()
        predecessor[e] = 0

        matches = np.all(data == predecessor[None, :], axis=1)

        if not matches.any():
            continue

        predecessor_index = int(np.argmax(matches))

        delta_MHN.append(
            comp(Ps_MHN[:, sample_index], Ps_MHN[:, predecessor_index])
        )
        delta_CBN.append(
            comp(Ps_CBN[:, sample_index], Ps_CBN[:, predecessor_index])
        )
        delta_baserate.append(
            comp(
                Ps_baserate[:, sample_index], Ps_baserate[:, predecessor_index]
            )
        )

delta_MHN = np.array(delta_MHN)
delta_CBN = np.array(delta_CBN)
delta_baserate = np.array(delta_baserate)

# print mean and 95% confidence interval around the mean
print(
    f"MHN: {np.mean(delta_MHN)} +- {1.96*np.std(delta_MHN, ddof=1)/np.sqrt(len(delta_MHN))}"
)
print(
    f"CBN: {np.mean(delta_CBN)} +- {1.96*np.std(delta_CBN, ddof=1)/np.sqrt(len(delta_CBN))}"
)
print(
    f"baserate: {np.mean(delta_baserate)} +- {1.96*np.std(delta_baserate, ddof=1)/np.sqrt(len(delta_baserate))}"
)
