# De-Novo Classification of Pan-Cancer Samples using large approximate Mutual Hazard Networks

This repository contains the data and information necessary to reproduce the results from the following paper:

De-Novo Classification of Pan-Cancer Samples using large approximate Mutual Hazard Networks  
Simon Pfahler, Andreas Lösch, Y. Linda Hu, Rudolf Schill, Stefan Hansch, Rainer Spang, Tilo Wettig

This Readme explains what the different parts of the repository are for and how they can be used.

## Dataset
The folder `dataset` provides the stratified dataset that was used for this analysis.
The data is from the MSK-CHORD dataset, and was accessed via cbioportal.

- `gene_panel.txt`: Comma-separated list of all genetic events in the dataset.
- `KM_Data.txt`: Survival data for all samples in the dataset.
- `sample_matrices`: Folder, contains the binary sample matrices for all cancer types individually, and the combined pan-cancer sample matrix that contains all cancer types.
- `test_data.txt`: Folder, small toy dataset that contains only the first 100 patients and the first 10 events in the dataset. To be used for testing only.

## Classification Results
The folder `classification_results` contains the classification obtained by the oMHN, CBN and base rate model for 13 groups. Each file contains the membership probabilities for all samples of the dataset in a space-separated matrix where each column represents a sample and each row represents a group.
Samples are sorted as in the files in `dataset`.

- `sample_Ps_oMHN_13groups.dat`: Membership probabilities obtained via the oMHN model
- `sample_Ps_CBN_13groups.dat`: Membership probabilities obtained via the CBN model
- `sample_Ps_baserate_13groups.dat`: Membership probabilities obtained via the base-rate model
- `classification_oMHN_13groups.dat`: Classification data obtained via the oMHN model
- `classification_CBN_13groups.dat`: Classification data obtained via the CBN model
- `classification_baserate_13groups.dat`: Classification data obtained via the base-rate model
- `thetas_oMHN`: Folder containing theta matrices for all 13 groups. Events are ordered as in `dataset/gene_panel.txt`, with observation rates in the last row.

## Validation of approximate MHN inference
The folder `validation` provides the code to validate the approximate MHN inference.

To generate the data shown in Fig. 1 and 2, run:
```
python generate_data.py
python test_contributions.py
python test_learning.py
python test_learning_score.py
```

- `generate_data.py`: Generate a dataset using the parameters from `parameters.py`. Data will be stored in `validation_data`.
- `parameters.py`: Set the total number of events, the number of samples, the maximum number of active events in the dataset, the mean baserate, the values of `dmax`, and the number of runs.
- `test_contributions.py`: Calculate score and gradient contributions for all thetas created. Results are stored in `logpTheta_entries` and `gradient_entries`.
- `test_learning.py`: Train exact and approximate MHNs for all datasets created. Results are stored in `learned_thetas`.
- `test_learning_score.py`: Calculate scores for all MHNs learned with `test_learning.py`.

## Progression
The folder `progression` contains a python script `progression.py` to generate figures for the evolution of the membership probabilities along a chain of events, like Fig. 4.

## oMHN classification
The folder `oMHN_classification` includes files to perform the classification into individual subgroups.
This code is written in python, and heavily uses the python library `fastmhn`, which can be installed [here](https://phygit.ur.de/physics/mhn/fastmhn).

As the classification takes a long time, it is strongly recommended run it on a compute cluster configured via [slurm](https://slurm.schedmd.com/overview.html).
If running using slurm, set `slurm=True` in `dataset/parameters.py`.
The code can run without slurm, but results on any meaningful dataset will take a long time (i.e. weeks).

To perform a classification, first update the values in `parameters.py` to the desired values, configure the `srun_*` slurm scripts (if applicable) and run `bash run_classification`.
The script runs in the background.

The overall progress is available in `log.json`, and the progress of the current iteration is available in `saved_aggregate.dat`.

After a successful run, classifications, theta matrices and aggregate information about every iteration are available in `sample_Ps`, `theta_it*` and `aggregates` respectively.

The following list briefly describes what all of the files are doing.

- `aggregate_CV_runs.py`: Check if cross-validation has found the optimal regularization strength. If not, it starts the next cross-validation runs.
- `cleanup.sh`: Remove all files that were created during a classification process.
- `get_models.py`: Obtain the final oMHN models of an EM iteration.
- `learn_approx_omhn_crossvalidated.py`: Obtain a cross-validated score for one group and regularization strength.
- `learn_approx_omhn.py`: Obtain an oMHN model for one group and regularization strength.
- `parameters.py`: Set the number of groups, the approximation parameter `dmax`, the number of EM iterations, the path to the dataset and if slurm is used.
- `run_classification`: Start a classification run.
- `split_samples.py`: Script to split the dataset into groups based on the current oMHN models.
- `srun_crossvalidation`: Run `learn_approx_omhn_crossvalidated.py` on a slurm cluster.
- `srun_model`: Run `learn_approx_omhn.py` on a slurm cluster.
- `srun_split_samples`: Run `split_samples.py` on a slurm cluster.
- `warden.py`: Train a classification. This script performs one E and M step, then restarts itself.

## CBN classification
The folder `CBN_classification` contains an R script `CBN_classification.r` that performs the classification based on the CBN method, given a number of clusters.

## base rate classification
The folder `baserate_classification` contains a python script `baserate_classification.py` that performs the classification based on the base rate method, for all cluster sizes from 2 to 50.
