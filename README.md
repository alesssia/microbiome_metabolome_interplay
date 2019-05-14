# README

These scripts are part of the manuscript: 

> Visconti A, *et al.*, *"Interplay between the human gut microbiome and host metabolism"*, bioXviv (2019), [https://www.biorxiv.org/content/10.1101/561787v1](https://www.biorxiv.org/content/10.1101/561787v1)

They describe the two custom analyses devised to test whether the microbiome is involved in the dialogue between the faecal and host systemic metabolism.

The scripts use age at the sample collection and sex as confounders, but any could be included by modifying the source code.

## Usage

These scripts are those used for the dialogue involving microbial metabolic pathways. The same scripts could be used to evaluate the dialogue involving species by simply loading the corresponding species-related files (the structure of the species/pathways file is the same, and described below). Example data are provided in the [`data`](./data) folder and described below.

Both scripts should be executed from the folder containing the data.

## Input data

- microbiome relative abundances and confounders: tab-separated file, describing species/metabolic pathways relative abundances. Abundances should have been previously arcsine square-root transformed, filtered for outliers using the Grubbs outlier test (significance threshold P=0.05), and standardised to have zero mean and unit variance. This file should also includes the confounders (age, sex)
  Columns should be:
  ```
  FID: family ID
  IID: individual ID
  age: age at sample collection
  sex: sex
  covariate1: other covariate
  ...
  covariateN: other covariate
  PWY.1: relative abundances of the first microbial metabolites pathway
  ...
  PWY.N: relative abundances of the Nth microbial metabolites pathway
  ```
- Blood metabolites: comma-separated log-transformed blood metabolite levels.
  Columns should be:
  ```
  FID: family ID
  IID: individual ID
  B1: log-transformed level of the first metabolite
  ...
  BN: log-transformed level of the Nth metabolite
  ```
- Faecal metabolites: comma-separated log-transformed faecal metabolite levels
  Columns should be:
  ```
  FID: family ID
  IID: individual ID
  F1: log-transformed level of the first metabolite
  ...
  FN: log-transformed level of the Nth metabolite
  ```
- tab-separated file describing the summary statistics of the association study between microbial pathways/species and faecal metabolites
  Columns should be:
  ```
  Pathway: associated pathway
  Faecal_metabolite: associated metabolite
  N: number of observation
  Beta: effect size
  SE: standard error
  P: P-value
  adj_P: adjusted P value
  ```
- tab-separated describing the summary statistics of the association study between microbial pathways/species and blood metabolites
  Columns should be:
  ```
  Pathway: associated pathway
  Blood_metabolite: associated metabolite
  N: number of observation
  Beta: effect size
  SE: standard error
  P: P-value
  adj_P: adjusted P value
  ```

## `simulation.R`

This script is used to test, for each pair of co-associated metabolites, whether a pathway is involved in the dialogue between faecal and blood metabolites. Indeed, in this case, these metabolites would be expected to be more strongly correlated in the presence of the pathway than in its absence. We used the missingness observed in our WMGS data to test this hypothesis via simulations.

To this aim, the script:
- selects all pairs of co-associated metabolites interacting with pathways with at least 30 missing observations
- builds 1,000 random datasets which included 1,000 pairs of metabolites matched by correlation and sample size to the original set of co-associated metabolites
- combines these new pairs with pathways having the same missingness pattern of the actually associated pathways
- uses these simulated datasets to assess the probability of observing an increased correlation between metabolites when the pathway is present in the co-associated metabolites compared to the matched pairs. 

Please note, that when using the example data, due to its small size, the matching is not possible, and the script does not terminate.

## `pgain.R`

This script is used to evaluate, for each pair of co-associated metabolites, its P-gain statistic, which allows determining whether the ratio between the two metabolites is more informative than the single metabolites alone, therefore suggesting the presence of a relationship between them [[DOI: 10.1186/1471-2105-13-120](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-120)]. 

To this aim, the script evaluates:
-  all the log ratios between each pair of co-associated metabolites (observed in at least 100 individuals with metagenomic data)
- the association with the single metabolites and with the obtained ratio with the specific pathway by fitting a linear mixed effect model in R, including age and sex as fixed effects, and family structure as a random effect. 
- the P-gain statistic as the ratio between the minimum P value obtained using the single metabolites alone and the P value obtained using their ratio

It has been previously observed that the magnitude of the P-gain statistic can be reduced by the increasing correlation between the metabolites and their ratio and increased by increasing sample size [[DOI: 10.1186/1471-2105-13-120](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-120)], two parameters which varied greatly in our dataset. 
Therefore, the script is also used to estimate a null distribution empirically using a conservative assumption of no interplay between metabolites associated with different pathways. 

To this aim, the script
- builds a null distribution of P-gain statistics using 100,000 pairs of randomly selected metabolites which were associated at 5% FDR with two different pathways by matching them 1-to-1 by correlation and sample size to the co-associated metabolite pairs. 

Please note, that when using the example data, due to its small size, the matching is not possible, and the script does not terminate.


## Example data

Example data are provided for microbial metabolic pathways.

Files are:
- community_pathways.tsv: microbial metabolic pathways relative abundances previously arcsine square-root transformed, filtered for outliers using the Grubbs outlier test (significance threshold P=0.05), and standardised to have zero mean and unit variance; this file includes also the covariates
- Blood_metabolites.csv: log-transformed blood metabolite levels
- Faecal_metabolites.csv: log-transformed faecal metabolite levels
- FaecalBloodCorrectAgeSexFam.tsv: faecal and blood metabolites corrected for confounders (age and sex)
- pathways_vs_blood.tsv: summary statistics of the association study between pathways and blood metabolites
- pathways_vs_faecal.tsv: summary statistics of the association study between pathways and faecal metabolites

Please note that, for the sake of the example, adjusted P values (adj_P) have been simulated.

