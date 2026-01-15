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

# >>> permute groups
mbs = data.sum(axis=1)


def sort_groups(Ps):
    mbs_groups = np.zeros(Ps.shape[0])
    nr_samples_groups = np.zeros(Ps.shape[0])
    for i in range(mbs.shape[0]):
        group = np.argmax(Ps[:, i])
        mbs_groups[group] += mbs[i]
        nr_samples_groups[group] += 1
    sorted_indices = np.argsort(mbs_groups / nr_samples_groups)[::-1]

    return Ps[sorted_indices]


Ps_MHN = sort_groups(Ps_MHN)
Ps_CBN = sort_groups(Ps_CBN)
Ps_baserate = sort_groups(Ps_baserate)
# <<< permute groups

delta_MHN = 0
delta_CBN = 0
delta_baserate = 0
nr_pairings = 0

# [((sample_a), (sample_b))]
steps_considered = list()


def comp(p1, p2):
    return np.linalg.norm(p1 - p2)


for sample_index in tqdm(range(Ps_MHN.shape[1])):
    for e in range(d):
        test_sample = data[sample_index]
        if test_sample[e] == 0:
            continue
        test_sample[e] = 0
        matches = np.all(data == test_sample[None, :], axis=1)

        if not matches.any():
            continue

        other_sample_index = int(np.argmax(matches))

        delta_MHN += comp(
            Ps_MHN[:, sample_index], Ps_MHN[:, other_sample_index]
        )
        delta_CBN += comp(
            Ps_CBN[:, sample_index], Ps_CBN[:, other_sample_index]
        )
        delta_baserate += comp(
            Ps_baserate[:, sample_index], Ps_baserate[:, other_sample_index]
        )
        nr_pairings += 1

        steps_considered.append((test_sample, data[sample_index]))

print(f"MHN: {delta_MHN / nr_pairings}")
print(f"CBN: {delta_CBN / nr_pairings}")
print(f"baserate: {delta_baserate / nr_pairings}")
