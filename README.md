# Discovering multiple antibiotic resistance phenotypes using diverse top-k subgroup list discovery

This repository contains all the needed resources in order to reproduce the research from the paper "Discovering multiple antibiotic resistance phenotypes using diverse top-k subgroup list discovery". For this, it is necessary to execute the files in the following order:

## 1. Initial query

The first step consists of executing a SQL query over the MIMIC-III database, which was previously stored in a PostgreSQL database. The [Initial-query/query.txt](Initial-query/query.txt) file contains this query, which produces a dataset in CSV format.

## 2. Preprocessing

After executing the SQL query, preprocessing is applied over the obtained dataset in order to generate the mining view. This preprocessing is included and detailed in the [Preprocessing/notebook.ipynb](Preprocessing/notebook.ipynb) file.

## 3. Experiments

Once having the mining view, it is necessary to execute the [Experiments/notebook.ipynb](Experiments/notebook.ipynb) file in order to prepare it for the experiments.

Next, the [Experiments/subgroups/mine_subgroups.py](Experiments/subgroups/mine_subgroups.py) file contains all the needed code in order to execute the [VLSD algorithm](https://doi.org/10.3390/a16060274) and mine the candidate subgroups.

After that, the format of the file generated above must be transformed to be able to execute our approach. This transformation is carried out with the [Experiments/phenotypes/transform_subgroups.py](Experiments/phenotypes/transform_subgroups.py) file.

Now, we can use this file with the candidate subgroups in order to mine diverse top-k phenotypes. This mining process is included in the [Experiments/phenotypes/mine_phenotypes.py](Experiments/phenotypes/mine_phenotypes.py) file.

When we have mined all the diverse top-k phenotypes, we select the best model using the code included in the [Experiments/phenotypes/select_best_model.ipynb](Experiments/phenotypes/select_best_model.ipynb) file.

Finally, the code to carry out the statistical validation of the best diverse top-k phenotypes is included in the [Experiments/statistical_validation/notebook.ipynb](Experiments/statistical_validation/notebook.ipynb) file.
