# Run SaintExpress from R on Windows or Linux

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
