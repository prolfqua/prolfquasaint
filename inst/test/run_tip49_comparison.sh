#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"

fixture_dir="${repo_root}/inst/SaintExpress/TIP49"
reference_bin="SAINTexpress-spc"
build_root="${SAINTEXPRESS_BUILD_ROOT:-${TMPDIR:-/tmp}/prolfquasaint-saintexpress-builds}"
build_301_dir="${build_root}/SAINTexpress-v3.0.1"
build_363_dir="${build_root}/SAINTexpress-v3.6.3"
native_301_bin="${build_301_dir}/SAINTexpress-spc"
native_363_bin="${build_363_dir}/SAINTexpress-spc"

build_native() {
  local source_dir="$1"
  local build_dir="$2"

  cmake -S "${source_dir}" -B "${build_dir}" \
    -DSAINTEXPRESS_TEST_DATA_DIR="${fixture_dir}" >/dev/null
  cmake --build "${build_dir}" >/dev/null
}

build_native "${repo_root}/inst/SAINTexpress-v3.0.1" "${build_301_dir}"
build_native "${repo_root}/inst/SAINTexpress-v3.6.3" "${build_363_dir}"

run_comparison() {
  local label="$1"
  local native_bin="$2"
  local native_dir="$3"
  local reference_dir="$4"
  local diff_file="$5"
  local report_file="$6"

  rm -rf "${native_dir}" "${reference_dir}" "${diff_file}" "${report_file}"
  mkdir -p "${native_dir}" "${reference_dir}"

  for input in inter.dat prey.dat bait.dat; do
    cp "${fixture_dir}/${input}" "${native_dir}/"
    cp "${fixture_dir}/${input}" "${reference_dir}/"
  done

  (
    cd "${native_dir}"
    "${native_bin}" inter.dat prey.dat bait.dat
  )

  (
    cd "${reference_dir}"
    docker run --rm --platform linux/amd64 \
      -v "${reference_dir}:/work" \
      -w /work \
      saintexpress:latest \
      "${reference_bin}" inter.dat prey.dat bait.dat
  )

  diff -u \
    "${reference_dir}/list.txt" \
    "${native_dir}/list.txt" \
    > "${diff_file}" || true

  if command -v quarto >/dev/null 2>&1; then
    (
      cd "${script_dir}"
      SAINTEXPRESS_REPORT_LABEL="${label}" \
      SAINTEXPRESS_NATIVE_LIST="$(basename "${native_dir}")/list.txt" \
      SAINTEXPRESS_REFERENCE_LIST="$(basename "${reference_dir}")/list.txt" \
        quarto render tip49_comparison.qmd \
          --output "$(basename "${report_file}")" \
          --execute-dir "${script_dir}"
    )
  else
    echo "quarto not found; skipping report render for ${label}" >&2
  fi

  echo "${label}"
  echo "Native output: ${native_dir}/list.txt"
  echo "Reference output: ${reference_dir}/list.txt"
  echo "Diff: ${diff_file}"
  echo "Report: ${report_file}"
}

run_comparison \
  "SAINTexpress 3.0.1 native spc vs packaged Linux spc" \
  "${native_301_bin}" \
  "${script_dir}/saintexpress-tip49-native" \
  "${script_dir}/saintexpress-tip49-reference" \
  "${script_dir}/saintexpress-tip49-reference-vs-native.diff" \
  "${script_dir}/tip49_comparison.html"

run_comparison \
  "SAINTexpress 3.6.3 native spc vs packaged Linux spc" \
  "${native_363_bin}" \
  "${script_dir}/saintexpress-363-tip49-native-spc" \
  "${script_dir}/saintexpress-363-tip49-reference-spc" \
  "${script_dir}/saintexpress-363-tip49-spc-reference-vs-native.diff" \
  "${script_dir}/tip49_comparison_363_spc.html"
