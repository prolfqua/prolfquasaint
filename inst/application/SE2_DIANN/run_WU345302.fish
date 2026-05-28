#!/usr/bin/env fish

set -l script (Rscript --vanilla -e "cat(system.file('application/bin/prolfqua_dea.sh', package = 'prolfquapp'))")
set -l base_dir (dirname (status --current-filename))/WU345302
set -l out_dir "$base_dir/prolfquapp_saint_run_tabs"

if test ! -x "$script"
    echo "prolfqua_dea.sh is not executable or was not found: $script" >&2
    exit 1
end

bash "$script" \
    --indir "$base_dir/out-DIANN" \
    --dataset "$base_dir/dataset_saint.csv" \
    --yaml "$base_dir/config.yaml" \
    --workunit 345302 \
    --software prolfquapp.DIANN \
    --outdir "$out_dir"
