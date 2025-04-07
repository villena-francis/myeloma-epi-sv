[![Integration test](https://github.com/cnio-bu/myeloma-epi-sv/actions/workflows/integration-test.yaml/badge.svg)](https://github.com/cnio-bu/myeloma-epi-sv/actions/workflows/integration-test.yaml)

# Myeloma epigenetics and SV analysis

## Running the pipeline

Due to dependencies on docker containers to run certain rules,
the Snakemake command needs to bind certain directories, e.g.:

```bash
snakemake \
    --sdm conda apptainer \
    --apptainer-args "-B /results_in_host/:/results/,/clairs_modes_in_host/:/models/,/logs_in_host/:/logs/" \
    --executor slurm \
    --conda-frontend conda # mamba fails when using together with apptainer
    ...
```
