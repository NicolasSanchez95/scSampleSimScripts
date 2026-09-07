# Simulation and all-hypothesis metrics (code availability)

These are companion scripts that generate the data used for the paper. Major analysis function live in the R package **postselect** ([epurdom/postselectPaper](https://github.com/epurdom/postselectPaper)). These scripts **install and load git tag** `tested-backup-20260904` (commit `ba61526`). They then make calls using SLURM to this analysis package to simulate scRNA-seq data with known DE structure, extract per-simulation master tables, and compute cluster-level and  experiment level metrics (FDP / power / imbalance).

Analysis functions live 

## Software

- R 4.x and Bioconductor packages used by `postselect` (including Seurat, muscat, harmony, batchelor, edgeR, scuttle, variancePartition)
- `remotes` and `BiocManager` (the first run installs tag `tested-backup-20260904` and any missing Bioconductor Imports such as `batchelor`)
- SLURM for the provided wrappers (`run_param_creator.sh`, `run_sims_extraction.sh`). Set `SLURM_PARTITION` and/or `SLURM_ACCOUNT` if your site requires them. There is no default partition or node list.

Optional: `export R_SCRIPT=/path/to/Rscript` if `Rscript` is not on `PATH`. Do not set `R_LIBS_USER` unless you want a non-default library.

## Input data

The SCE object is not pushed. It is assumed to be placed at:

```text
simulation_scripts/data/filtered_sce_data.Rda
```

Required `colData` columns: `group`, `sample`, `batch`, `cell_type`. Assay `logcounts` is used during simulation setup.

## Run order

From the repository root (or with `--repo-root` pointing here).

**1. Parameter files and 50 simulations** (one DE-gene count, default 100).

The param creator is a heavy R job (Harmony / Leiden on the SCE). Submit it with SLURM rather than running it on a login node. Flags after the script name are passed through:

```bash
cd simulation_scripts
mkdir -p logs
sbatch -p PARTITION run_param_creator.sh --sim_prefix augSimLeiden \
  --sim_type augData --cluster_type leiden --n_sims 50 --num_de_genes 100
```

That job writes param files, then submits the 50-simulation array on the same partition. Logs: `simulation_scripts/logs/param_creator_<jobid>.out`.

For the k-means hierarchy example:

```bash
sbatch -p PARTITION run_param_creator.sh --sim_prefix augSimHierarchy \
  --sim_type augData --cluster_type kmeans --n_sims 50 --num_de_genes 100
```

`param_creator_augSimLeiden.R` and `param_creator_augSimHierarchy.R` are the two example parameter settings (Leiden vs k-means hierarchy `cut_at` / LFC) while sim_analysis_pipeline.R does the simulation and post-simulation analysis with the different phenotype removals. 

**2. Master tables** (after simulations finish):

By doing a regex search in the results folder for the provided id-pattern, this will create a cohort (list of ids) and generate a master table with all information about all tested and simulated gene-cluster hypotheses. These are large files.

```bash
cd master_table_scripts
./run_sims_extraction.sh --cohort-name demo50 --id-pattern "augSimLeiden_"
```

**3. Cluster and Experiment metrics:**

```bash
cd decision_parameters
./run_allHyp_metrics.sh --cohort-name demo50
```

The first R job that needs `postselect` installs tag `tested-backup-20260904` if the library SHA does not match `ba61526`. Later array tasks reuse that install.

## Outputs


| Stage         | Location                                                                                                    |
| ------------- | ----------------------------------------------------------------------------------------------------------- |
| Simulation    | `simulation_scripts/results/` (`anls_info/`, `de_final/`, `sims_info/`, `de_configs/`)                      |
| Master tables | `master_table_scripts/results/<cohort>/master_pwr_table_*.rds`                                              |
| Cluster level metrics  | `decision_parameters/results/<cohort>/conf_metrics_tested_*.rds` and `conf_metrics_tested_by_cluster_*.rds` |


Shared settings are in `config.sh` (repo root, `POSTSELECT_REF`, `POSTSELECT_SHA`). The same pin is in `R/ensure_postselect.R`.