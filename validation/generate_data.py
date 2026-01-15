import os

import fastmhn
import numpy as np

from parameters import *

np.random.seed(42)

os.makedirs("validation_data")

nr_theta = 0
while nr_theta < nr_thetas:

    # create ground truth theta
    theta = fastmhn.utility.generate_theta(
        d,
        base_rate_loc=mean_baserate,
        base_rate_scale=1,
        influence_loc=0,
        influence_scale=1,
        sparsity=0.9,
    )

    # save ground truth model
    np.savetxt(f"validation_data/theta{nr_theta}.dat", theta)

    # generate data from ground truth model
    data_raw = fastmhn.utility.generate_data(theta, 4 * nr_samples)

    # make sure data has no samples with >= 25 events, as that would be
    # problematic for the exact calculation
    data = np.zeros((nr_samples, d), dtype=np.int32)
    index_data = 0
    for index_data_raw, patient in enumerate(data_raw):
        if index_data >= nr_samples:
            print(
                f"Info: Skipped {index_data_raw-index_data}"
                f" samples with >={max_active_events} mutations."
            )
            break
        if np.sum(patient) <= max_active_events:
            data[index_data] = patient
            index_data += 1

    if index_data != nr_samples:
        print(
            f"Warning: Dataset {nr_theta} contains only"
            f" {index_data} samples!"
        )

    print(f"Dataset {nr_theta}")
    print(f"\tNr of samples: {data.shape[0]}")
    print(f"\tAvg mut burden: {np.mean(np.sum(data, axis=1))}")
    print(f"\tMax mut burden: {np.max(np.sum(data, axis=1))}")

    # check if the dataset is valid
    if np.any(np.sum(data, axis=0) == 0):
        print(f"Info: Dataset {nr_theta} is invalid, regenerating!")
    else:
        print(f"Info: Dataset {nr_theta} sucessfully generated.")
        # save dataset
        np.savetxt(f"validation_data/data{nr_theta}.dat", data)
        nr_theta += 1
