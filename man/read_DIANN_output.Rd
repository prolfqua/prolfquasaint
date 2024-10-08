% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_DIANN_output.R
\name{read_DIANN_output}
\alias{read_DIANN_output}
\title{DIANN helper functions}
\usage{
read_DIANN_output(
  diann.path,
  fasta.file,
  nrPeptides = 2,
  q_value = 0.01,
  isUniprot = TRUE,
  rev = "REV_"
)
}
\arguments{
\item{diann.path}{path to diann-output.tsv}

\item{fasta.file}{path to fasta file}

\item{nrPeptides}{peptide threshold}

\item{Q.Value}{q value threshold}
}
\description{
DIANN helper functions
}
\examples{
\dontrun{
x <- prolfquapp::get_DIANN_files("inst/application/DIANN/2517219/")
xd <- read_DIANN_output(x$data, x$fasta)
debug(read_DIANN_output)
}
}
