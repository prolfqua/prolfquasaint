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
    id_col = "protein_Id" ){

  plengths <- data.frame(id = names(fasta) , protein.length = vapply(fasta, stringr::str_length, integer(1)))
  byx <- "id"
  names(byx) <- id_col
  intdata <- dplyr::left_join(intdata, plengths , by = byx)
  intdata$protein.length[is.na(intdata$protein.length)] <- as.integer(mean(intdata$protein.length, na.rm = TRUE))
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
protein_2localSaint <- function(xx,
                                quantcolumn = "mq.protein.intensity",
                                proteinID = "protein_Id",
                                geneNames  = proteinID,
                                proteinLength = "protein.length",
                                IP_name = "raw.file",
                                baitCol = "bait",
                                CorTCol = "CorT"
){
  reqcolumns <- c(quantcolumn,proteinID,geneNames,proteinLength,IP_name,baitCol,CorTCol)
  if ( !all(reqcolumns %in% colnames(xx)) ) {

    stop("columns not found ", paste0(reqcolumns[which(!reqcolumns %in% colnames(xx))]))
  }
  res <- list()
  bait <- xx |> dplyr::select(!!!rlang::syms(c(IP_name,baitCol,CorTCol)))
  bait <- dplyr::distinct(bait)
  res$bait <- bait
  prey <- xx |> dplyr::select(proteinID = !!rlang::sym(proteinID),
                              proteinLength = !!rlang::sym(proteinLength),
                              geneNames = !!rlang::sym(geneNames))
  prey <- dplyr::distinct(prey)
  res$prey <- prey
  inter <- xx |>
    dplyr::select(!!!rlang::syms(c(IP_name,
                            baitCol,
                            proteinID ,
                            quantcolumn))) |>
    dplyr::filter(!!rlang::sym(quantcolumn) > 0)
  res$inter <- inter
  res <- res[c("inter","prey","bait")]
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
ensure_saintexpress_docker_image <- function(image_name = "saintexpress:latest") {
  docker_path <- Sys.which("docker")
  if (docker_path == "") {
    stop(
      "Docker is not installed or not in PATH. ",
      "Install Docker Desktop from https://www.docker.com/products/docker-desktop/"
    )
  }

  check <- system2(
    "docker", args = c("image", "inspect", image_name),
    stdout = FALSE, stderr = FALSE
  )
  if (check == 0) {
    message("Docker image '", image_name, "' found.")
    return(invisible(TRUE))
  }

  message("Building Docker image '", image_name, "' (one-time setup)...")

  build_dir <- file.path(tempdir(), "saintexpress_build")
  dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)

  dockerfile_path <- prolfqua::find_package_file("prolfquasaint", "docker/Dockerfile")
  bin_dir <- prolfqua::find_package_file("prolfquasaint", "SaintExpress/bin/Linux64")
  file.copy(dockerfile_path, build_dir, overwrite = TRUE)
  file.copy(file.path(bin_dir, "SAINTexpress-spc"), build_dir, overwrite = TRUE)
  file.copy(file.path(bin_dir, "SAINTexpress-int"), build_dir, overwrite = TRUE)

  result <- system2(
    "docker",
    args = c(
      "build", "--platform", "linux/amd64",
      "-t", image_name, build_dir
    ),
    stdout = TRUE, stderr = TRUE
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


#' Run SAINTexpress analysis
#'
#' Executes SAINTexpress on prepared input data. On macOS, uses Docker with
#' Rosetta to run the Linux binary transparently (requires Docker Desktop).
#' On Linux and Windows, uses the native binary.
#'
#' @export
#' @rdname saintExpress
#' @param si output of protein_2localSaint function
#' @param filedir where to store the saint express inputs
#' @param spc if TRUE spectral counts, if FALSE intensities (see SAINTexpress documentation)
#' @param CLEANUP TRUE remove all files generated by SAINTexpress
#' @param use_docker logical or NULL. NULL (default) auto-detects: uses Docker
#'   on macOS, native binary elsewhere. TRUE forces Docker. FALSE forces native.
runSaint <- function(si,
                     filedir = getwd(),
                     spc = TRUE,
                     CLEANUP = TRUE,
                     use_docker = NULL) {
  stopifnot(names(si) == c("inter", "prey", "bait"))

  filedir <- normalizePath(filedir, mustWork = TRUE)

  paths <- character(3)
  for (i in seq_along(si)) {
    filen <- file.path(filedir, paste0(names(si)[i], ".txt"))
    paths[i] <- filen
    message(filen)
    readr::write_tsv(si[[i]], file = filen, col_names = FALSE)
  }

  sysname <- Sys.info()["sysname"]
  if (is.null(use_docker)) {
    use_docker <- (sysname == "Darwin")
  }

  binary_name <- if (spc) "SAINTexpress-spc" else "SAINTexpress-int"

  if (use_docker) {
    # --- Docker execution (macOS default, or forced) ---
    ensure_saintexpress_docker_image()
    container_dir <- "/data"
    container_paths <- file.path(container_dir, basename(paths))

    docker_args <- c(
      "run", "--rm",
      "--platform", "linux/amd64",
      "-v", paste0(filedir, ":", container_dir),
      "-w", container_dir,
      "saintexpress:latest",
      binary_name,
      container_paths
    )

    out <- system2(
      "docker", args = docker_args,
      stdout = TRUE, stderr = TRUE, wait = TRUE
    )
    listFile <- file.path(filedir, "list.txt")

  } else if (sysname == "Windows") {
    exeS2 <- prolfqua::find_package_file(
      "prolfquasaint",
      paste0("SaintExpress/bin/Windows64/", binary_name, ".exe")
    )
    out <- system2(
      exeS2, args = paths,
      stdout = TRUE, stderr = TRUE, wait = TRUE, minimized = TRUE
    )
    listFile <- file.path(getwd(), "list.txt")

  } else if (sysname == "Linux") {
    exeS2 <- prolfqua::find_package_file(
      "prolfquasaint",
      paste0("SaintExpress/bin/Linux64/", binary_name)
    )
    out <- system2(
      exeS2, args = paths,
      stdout = TRUE, stderr = TRUE, wait = TRUE
    )
    listFile <- file.path(getwd(), "list.txt")

  } else {
    stop(
      "System ", sysname, " not supported natively. ",
      "Set use_docker = TRUE to use Docker."
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

