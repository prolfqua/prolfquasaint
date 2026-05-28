# TODO: Add WU345302 SAINTexpress prolfquapp Fish Runner

## Goal

Add a small fish script beside the `WU345302` folder that records the `prolfqua_dea.sh` command and parameters used to run the prolfquapp SAINT workflow.

## Plan

1. Use the existing prolfquapp-compatible `WU345302/config.yaml` and `dataset_saint.csv`.
2. Create `inst/application/SE2_DIANN/run_WU345302.fish` as a sibling of the workunit folder.
3. Keep output directed to `WU345302/prolfquapp_saint_run_tabs`.
4. Make the script executable and verify its contents.
