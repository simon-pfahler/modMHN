import matplotlib.pyplot as plt
import numpy as np

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

active_events = (
    input("CSV list of active event names, leave empty for an example: ")
    .strip()
    .split(",")
)
if active_events == [""]:
    active_events = (
        "PIK3CA,TP53,APC,KRAS,SOX9,FBXW7,TGFBR2,MAP3K1".strip().split(",")
    )

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

Ps_baserate_progression = np.zeros((groups, len(active_events) + 1))
Ps_MHN_progression = np.zeros((groups, len(active_events) + 1))
Ps_CBN_progression = np.zeros((groups, len(active_events) + 1))


def find_sample(active_events):
    global data_progression

    sample = np.zeros(d)
    for e in active_events:
        sample[eventnames.index(e)] = 1

    index = np.where(np.all(data == sample, axis=1))[0]
    if index.size:
        print(f"Sample found, first at index {index[0]}")
        Ps_baserate_progression[:, len(active_events)] = Ps_baserate[
            :, index[0]
        ]
        Ps_MHN_progression[:, len(active_events)] = Ps_MHN[:, index[0]]
        Ps_CBN_progression[:, len(active_events)] = Ps_CBN[:, index[0]]
    else:
        print(f"Sample not found!")
        Ps_baserate_progression[:, len(active_events)] = np.nan
        Ps_MHN_progression[:, len(active_events)] = np.nan
        Ps_CBN_progression[:, len(active_events)] = np.nan
    return sample


sample = None
for i in range(len(active_events) + 1):
    sample = find_sample(active_events[:i])

order_MHN = np.argsort(Ps_MHN_progression.sum(axis=1))
Ps_MHN_progression = Ps_MHN_progression[order_MHN]
group_names_MHN = [group_names_MHN[i] for i in order_MHN]

order_CBN = np.argsort(Ps_CBN_progression.sum(axis=1))
Ps_CBN_progression = Ps_CBN_progression[order_CBN]
group_names_CBN = [group_names_CBN[i] for i in order_CBN]

order_baserate = np.argsort(Ps_baserate_progression.sum(axis=1))
Ps_baserate_progression = Ps_baserate_progression[order_baserate]
group_names_baserate = [group_names_baserate[i] for i in order_baserate]

x = np.arange(Ps_baserate_progression.shape[1])

fig, ax = plt.subplots(1, 3)

ax[0].stackplot(x, Ps_MHN_progression, labels=group_names_MHN)
ax[0].set_title("oMHN")
ax[0].set_ylim(0, 1)
ax[0].set_xlim(0, len(active_events))
ax[0].set_xticks(
    list(range(len(active_events) + 1)), ["none"] + active_events, rotation=90
)
ax[0].legend(loc="lower left")

ax[1].stackplot(x, Ps_CBN_progression, labels=group_names_CBN)
ax[1].set_title("CBN")
ax[1].set_ylim(0, 1)
ax[1].set_xlim(0, len(active_events))
ax[1].set_xticks(
    list(range(len(active_events) + 1)), ["none"] + active_events, rotation=90
)
ax[1].legend(loc="lower left")

ax[2].stackplot(x, Ps_baserate_progression, labels=group_names_baserate)
ax[2].set_title("base-rate")
ax[2].set_ylim(0, 1)
ax[2].set_xlim(0, len(active_events))
ax[2].set_xticks(
    list(range(len(active_events) + 1)), ["none"] + active_events, rotation=90
)
ax[2].legend(loc="lower left")

plt.tight_layout()
plt.show()


def comp(p1, p2):
    return np.sum(p1 * np.log(p1 / p2))  # KL divergence


for i in range(len(active_events)):
    print(
        f"{active_events[i]}: {comp(Ps_MHN_progression[:,i+1],Ps_MHN_progression[:,i])}, {comp(Ps_CBN_progression[:,i+1],Ps_CBN_progression[:,i])}, {comp(Ps_baserate_progression[:,i+1],Ps_baserate_progression[:,i])}"
    )
