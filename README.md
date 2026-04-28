# Run SaintExpress from R on Windows or Linux

This package wraps SAINTexpress for protein interaction scoring. The packaged
SAINTexpress 3.6.3 source archive was downloaded from the `saint-apms`
SourceForge project:

- SourceForge files page: https://sourceforge.net/projects/saint-apms/files/
- Archive: `SAINTexpress_v3.6.3__2018-03-09.tar.gz`
- SourceForge listing date: 2018-03-15

The extracted and cleaned C++ source is maintained separately under
`inst/SAINTexpress-v3.6.3`. It builds `SAINTexpress-spc` and
`SAINTexpress-int` with CMake against system Boost.Program_options and NLopt,
without vendored Boost/NLopt sources or precompiled SAINTexpress binaries.

On macOS:

```bash
brew install boost nlopt
cmake -S inst/SAINTexpress-v3.6.3 -B inst/SAINTexpress-v3.6.3/build
cmake --build inst/SAINTexpress-v3.6.3/build
ctest --test-dir inst/SAINTexpress-v3.6.3/build --output-on-failure
```

The native build has been compared against the packaged Linux SAINTexpress
3.6.3 binary using the TIP49 fixture in `inst/test/run_tip49_comparison.sh`.
The comparison is very close but not byte-identical, with small numerical
differences in a subset of scores.

How to install.

```bash

export R_LIBS_SITE="/scratch/PROLFQUA/r-site-library/"
R --vanilla << EOF
.libPaths()
install.packages(c("remotes","seqinr", "prozor","logger"), repos = "https://stat.ethz.ch/CRAN/")
remotes::install_gitlab("wolski/prolfquadata", host="gitlab.bfabric.org")
remotes::install_github("https://github.com/fgcz/prolfqua", build_vignettes = TRUE, dependencies = TRUE)
remotes::install_github("https://github.com/prolfqua/prolfquapp", dependencies = TRUE)
remotes::install_github("https://github.com/prolfqua/prolfquasaint", dependencies = TRUE)
EOF

```
