#' Add protein lengths from fasta file to data frame (id_col - protein id column.)
#' @rdname saintExpress
#' @param intdata data.frame
#' @param fasta list of sequences created with \code{\link[seqinr]{read.fasta}}
#' @param id_col column with protein ids/accessions.
#' @export
#'
add_protein_lengths <- function(
  intdata,
  fasta,
  id_col = "protein_Id"
) {
  plengths <- data.frame(
    id = names(fasta),
    protein.length = vapply(fasta, stringr::str_length, integer(1))
  )
  byx <- "id"
  names(byx) <- id_col
  intdata <- dplyr::left_join(intdata, plengths, by = byx)
  intdata$protein.length[is.na(intdata$protein.length)] <- as.integer(mean(
    intdata$protein.length,
    na.rm = TRUE
  ))
  return(intdata)
}

#' Convert tidy table with protein quants into SaintExpress compatible inputs
#' @param xx data.frame in long format
#' @export
#' @rdname saintExpress
#' @param quantcolumn intensity column
#' @param proteinID protein accession
#' @param geneNames column with gene names
#' @param proteinLength column with protein lengths
#' @param IP_name raw.file
#' @param baitCol column with bait definition (condition)
#' @param CorTCol is it control or TRUE (SaintExpress speach)
#' @examples
#'
#' bb <- prolfqua::prolfqua_data('data_IonstarProtein_subsetNorm')
#' bb$config <- bb$config$clone(deep = TRUE)
#' xx <- prolfqua::LFQData$new(bb$data, bb$config)
#' exampleDat <- xx$data |> dplyr::mutate(CorT = dplyr::case_when(dilution. == "a" ~ "C", TRUE ~ "T"))
#' # sample protein lengths
#'
#' tmp <- data.frame(protein_Id = unique(exampleDat$protein_Id))
#' tmp$proteinLength <- as.integer(runif(nrow(tmp), min = 150, max = 2500))
#' exampleDat <- dplyr::inner_join(tmp, exampleDat)
#' #undebug(protein_2localSaint)
#' res <- protein_2localSaint(exampleDat,quantcolumn = "medpolish",
#'                    proteinID = "protein_Id",
#'                    proteinLength = "proteinLength",
#'                    IP_name = "raw.file",
#'                    baitCol = "dilution.",
#'                    CorTCol = "CorT"
#'                    )
#'
#' stopifnot(names(res) == c( "inter", "prey",  "bait"))
#' if (Sys.info()["sysname"] == "Darwin" && Sys.which("docker") == "") {
#'   testthat::expect_error(runSaint(res, filedir = tempdir()))
#' } else {
#'   data_SAINTe_output <- runSaint(res, filedir = tempdir())
#' }
#'
protein_2localSaint <- function(
  xx,
  quantcolumn = "mq.protein.intensity",
  proteinID = "protein_Id",
  geneNames = proteinID,
  proteinLength = "protein.length",
  IP_name = "raw.file",
  baitCol = "bait",
  CorTCol = "CorT"
) {
  reqcolumns <- c(
    quantcolumn,
    proteinID,
    geneNames,
    proteinLength,
    IP_name,
    baitCol,
    CorTCol
  )
  if (!all(reqcolumns %in% colnames(xx))) {
    stop(
      "columns not found ",
      paste0(reqcolumns[which(!reqcolumns %in% colnames(xx))])
    )
  }
  res <- list()
  bait <- xx |> dplyr::select(!!!rlang::syms(c(IP_name, baitCol, CorTCol)))
  bait <- dplyr::distinct(bait)
  res$bait <- bait
  prey <- xx |>
    dplyr::select(
      proteinID = !!rlang::sym(proteinID),
      proteinLength = !!rlang::sym(proteinLength),
      geneNames = !!rlang::sym(geneNames)
    )
  prey <- dplyr::distinct(prey)
  res$prey <- prey
  inter <- xx |>
    dplyr::select(
      !!!rlang::syms(c(IP_name, baitCol, proteinID, quantcolumn))
    ) |>
    dplyr::filter(!!rlang::sym(quantcolumn) > 0)
  res$inter <- inter
  res <- res[c("inter", "prey", "bait")]
  return(res)
}


#' Ensure the SAINTexpress Docker image exists, building it if necessary
#'
#' Called automatically by \code{\link{runSaint}} on macOS. Checks for a local
#' Docker image and builds it from the shipped Dockerfile + Linux binaries
#' on first use (~10 seconds one-time cost).
#'
#' @param image_name Docker image name (default \code{"saintexpress:latest"})
#' @return \code{TRUE} invisibly on success
#' @keywords internal
ensure_saintexpress_docker_image <- function(
  image_name = "saintexpress:latest"
) {
  docker_path <- Sys.which("docker")
  if (docker_path == "") {
    stop(
      "Docker is not installed or not in PATH. ",
      "Install Docker Desktop from https://www.docker.com/products/docker-desktop/"
    )
  }

  check <- system2(
    "docker",
    args = c("image", "inspect", image_name),
    stdout = FALSE,
    stderr = FALSE
  )
  if (check == 0) {
    message("Docker image '", image_name, "' found.")
    return(invisible(TRUE))
  }

  message("Building Docker image '", image_name, "' (one-time setup)...")

  build_dir <- file.path(tempdir(), "saintexpress_build")
  dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)

  dockerfile_path <- prolfqua::find_package_file(
    "prolfquasaint",
    "docker/Dockerfile"
  )
  bin_dir <- prolfqua::find_package_file(
    "prolfquasaint",
    "SaintExpress/bin/Linux64"
  )
  file.copy(dockerfile_path, build_dir, overwrite = TRUE)
  file.copy(file.path(bin_dir, "SAINTexpress-spc"), build_dir, overwrite = TRUE)
  file.copy(file.path(bin_dir, "SAINTexpress-int"), build_dir, overwrite = TRUE)

  result <- system2(
    "docker",
    args = c(
      "build",
      "--platform",
      "linux/amd64",
      "-t",
      image_name,
      build_dir
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop(
      "Failed to build Docker image. Is Docker Desktop running?\n",
      paste(result, collapse = "\n")
    )
  }

  message("Docker image '", image_name, "' built successfully.")
  invisible(TRUE)
}

.saintexpress_package_file <- function(...) {
  rel <- file.path(...)
  installed <- system.file(rel, package = "prolfquasaint")
  if (nzchar(installed)) {
    return(installed)
  }
  source_tree <- file.path(getwd(), "inst", rel)
  if (file.exists(source_tree)) {
    return(normalizePath(source_tree, mustWork = TRUE))
  }
  parent_source_tree <- file.path(getwd(), "..", "inst", rel)
  if (file.exists(parent_source_tree)) {
    return(normalizePath(parent_source_tree, mustWork = TRUE))
  }
  ""
}

.saintexpress_executable <- function(
  binary_name,
  sysname = Sys.info()[["sysname"]]
) {
  candidates <- switch(
    sysname,
    "Darwin" = c(
      .saintexpress_package_file("SaintExpress", "bin", "Darwin", binary_name),
      .saintexpress_package_file("SaintExpress", "bin", "macOS", binary_name),
      .saintexpress_package_file("SAINTexpress-v3.6.3", "build", binary_name)
    ),
    "Windows" = c(
      .saintexpress_package_file(
        "SaintExpress",
        "bin",
        "Windows64",
        paste0(binary_name, ".exe")
      )
    ),
    "Linux" = c(
      .saintexpress_package_file("SaintExpress", "bin", "Linux64", binary_name)
    ),
    character()
  )
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[
    file.exists(candidates) & file.access(candidates, mode = 1) == 0
  ]
  if (!length(candidates)) {
    return("")
  }
  candidates[[1]]
}

.saintexpress_run_native <- function(
  exe,
  paths,
  filedir,
  sysname = Sys.info()[["sysname"]]
) {
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(filedir)
  args <- list(
    command = exe,
    args = paths,
    stdout = TRUE,
    stderr = TRUE,
    wait = TRUE
  )
  if (identical(sysname, "Windows")) {
    args$minimized <- TRUE
  }
  out <- do.call(system2, args)
  list(out = out, listFile = file.path(filedir, "list.txt"))
}


#' Run SAINTexpress analysis
#'
#' Executes SAINTexpress on prepared input data. Uses a native executable when
#' one is available for the current platform. On macOS, Docker can still be
#' forced with \code{use_docker = TRUE}.
#'
#' @export
#' @rdname saintExpress
#' @param si output of protein_2localSaint function
#' @param filedir where to store the saint express inputs
#' @param spc if TRUE spectral counts, if FALSE intensities (see SAINTexpress documentation)
#' @param CLEANUP TRUE remove all files generated by SAINTexpress
#' @param use_docker logical or NULL. NULL (default) uses a native executable
#'   when available and falls back to Docker on macOS. TRUE forces Docker.
#'   FALSE forces native execution.
#' @param engine execution engine. \code{"binary"} runs bundled SAINTexpress;
#'   \code{"r"} runs the R implementation.
#' @param optimizer optimizer for \code{engine = "r"}. \code{"base"} uses
#'   \code{\link[stats]{optim}}; \code{"nloptr"} uses NLopt COBYLA when
#'   \pkg{nloptr} is installed.
runSaint <- function(
  si,
  filedir = getwd(),
  spc = TRUE,
  CLEANUP = TRUE,
  use_docker = NULL,
  engine = c("binary", "r"),
  optimizer = c("base", "nloptr")
) {
  stopifnot(names(si) == c("inter", "prey", "bait"))
  engine <- match.arg(engine)
  optimizer <- match.arg(optimizer)

  if (engine == "r") {
    if (isTRUE(spc)) {
      res <- saint_spc_r(si, optimizer = optimizer)
      implementation <- "spectral-count"
    } else {
      res <- saint_int_r(si, optimizer = optimizer)
      implementation <- "intensity"
    }
    return(list(
      listFile = data.frame(listFile = NA_character_),
      list = res,
      out = data.frame(
        out = paste("R SAINTexpress", implementation, "implementation")
      )
    ))
  }

  filedir <- normalizePath(filedir, mustWork = TRUE)

  paths <- character(3)
  for (i in seq_along(si)) {
    filen <- file.path(filedir, paste0(names(si)[i], ".txt"))
    paths[i] <- filen
    message(filen)
    readr::write_tsv(si[[i]], file = filen, col_names = FALSE)
  }

  sysname <- Sys.info()[["sysname"]]
  binary_name <- if (spc) "SAINTexpress-spc" else "SAINTexpress-int"
  native_exe <- .saintexpress_executable(binary_name, sysname = sysname)
  if (is.null(use_docker)) {
    use_docker <- identical(sysname, "Darwin") && !nzchar(native_exe)
  }

  if (use_docker) {
    # --- Docker execution (macOS default, or forced) ---
    ensure_saintexpress_docker_image()
    container_dir <- "/data"
    container_paths <- file.path(container_dir, basename(paths))

    docker_args <- c(
      "run",
      "--rm",
      "--platform",
      "linux/amd64",
      "-v",
      paste0(filedir, ":", container_dir),
      "-w",
      container_dir,
      "saintexpress:latest",
      binary_name,
      container_paths
    )

    out <- system2(
      "docker",
      args = docker_args,
      stdout = TRUE,
      stderr = TRUE,
      wait = TRUE
    )
    listFile <- file.path(filedir, "list.txt")
  } else if (nzchar(native_exe)) {
    native <- .saintexpress_run_native(
      native_exe,
      paths = paths,
      filedir = filedir,
      sysname = sysname
    )
    out <- native$out
    listFile <- native$listFile
  } else {
    stop(
      "No native ",
      binary_name,
      " executable found for ",
      sysname,
      ". ",
      if (identical(sysname, "Darwin")) {
        "Build inst/SAINTexpress-v3.6.3 with CMake or set use_docker = TRUE."
      } else {
        "Install a platform binary or set use_docker = TRUE where Docker is supported."
      }
    )
  }

  message(cat(out, sep = "\n"))
  Sys.sleep(2)

  res <- utils::read.csv(file = listFile, sep = "\t")
  if (CLEANUP) {
    if (!file.remove(listFile)) {
      warning("can't remove ", listFile)
    }
    file.remove(paths)
  }
  res <- list(
    listFile = data.frame(listFile = listFile),
    list = res,
    out = data.frame(out = out)
  )
  return(res)
}
