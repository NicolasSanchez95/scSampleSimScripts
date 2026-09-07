# Master table extraction (cohorts + SLURM jobs)

This folder builds a **cohort** (simulation IDs) from analysis outputs, then
submits SLURM jobs that write `master_pwr_table_<phen_type>.rds`. 

Here <phen_type> will represent the used phenotype removal method.

## What to change

- `--cohort-name`: label used in filenames and output folders
- `--id-pattern`: `list.files()` regex over `anls_info/pa_de_int_*.rds`
- `PHENOTYPES` in `run_sims_extraction.sh`: phenotype-removal prefixes to extract

## Paths

- `--sim-results-dir`: simulation outputs; must contain `anls_info/` and `de_final/`
(default: `<repo-root>/simulation_scripts/results`)
- `--cohort-data-dir`: where `cohort_ids_<cohort-name>.rds` is written
(default: `<repo-root>/master_table_scripts/data`)
- `--results-base`: master-table outputs
(default: `<repo-root>/master_table_scripts/results`)



## Entrypoint

```bash
cd master_table_scripts
./run_sims_extraction.sh \
  --cohort-name demo50 \
  --id-pattern "augSimLeiden_"
```

Outputs:

- `<cohort-data-dir>/cohort_ids_<cohort-name>.rds`
- `<results-base>/<cohort-name>/master_pwr_table_<phen_type>.rds`

