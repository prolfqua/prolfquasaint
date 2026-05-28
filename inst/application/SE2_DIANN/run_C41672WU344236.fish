#!/usr/bin/env fish

set -l script (Rscript --vanilla -e "cat(system.file('application/bin/prolfqua_dea.sh', package = 'prolfquapp'))")
set -l base_dir (dirname (status --current-filename))/C41672WU344236
set -l out_dir "$base_dir/prolfquapp_saint_run_tabs"

if test ! -x "$script"
    echo "prolfqua_dea.sh is not executable or was not found: $script" >&2
    exit 1
end

bash "$script" \
    --indir "$base_dir" \
    --dataset "$base_dir/dataset.csv" \
    --yaml "$base_dir/config_prolfquapp_saint.yaml" \
    --workunit 344236 \
    --software prolfquapp.DIANN \
    --outdir "$out_dir"
