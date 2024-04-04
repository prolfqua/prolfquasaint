# Author : Witold Wolski <wew@fgcz.ethz.ch>
# compatible with prolfqua 3.0.0 release available from https://github.com/wolski/prolfqua/releases/tag/v0.2.9


prolfquasaint::copy_SAINT_express(run_script = FALSE)
# Read b-fabric related information
yml <- yaml::read_yaml("config.yaml")

BFABRIC <- prolfquasaint::get_params_Bfabric(yml)





# list with data used with the markdown report
REPORTDATA <- prolfquasaint::apparams_Bfabric(yml)
REPORTDATA$nrPeptides <- 2
ZIPDIR = paste0("C",BFABRIC$orderID,"WU", BFABRIC$workunitID, "NRpep", REPORTDATA$nrPeptides )
dir.create(ZIPDIR)

# Prefix for exported files
treat <- "DIANN_"
# load data

dataset <- dir(".", pattern = 'dataset.csv', full.names = TRUE, recursive = TRUE)[1]
annot <- read.csv(dataset)
annot <- data.frame(lapply(annot, as.character))
annot <- annot |> dplyr::mutate(
  raw.file = gsub("^x|.d.zip$|.raw$","",
                  (basename(annot$Relative.Path))
  )
)

annotation <- annot
colnames(annotation) <- tolower(make.names(colnames(annotation)))

path = "."
diann.path <- grep("report\\.tsv$|diann-output\\.tsv", dir(path = path, recursive = TRUE), value = TRUE)
fasta.files <- grep("*\\.fasta$|*\\.fas$", dir(path = path, recursive = TRUE), value = TRUE)

if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
  fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
}

if (length(fasta.files) == 0) {
  logger::log_error("No fasta file found!")
  stop()
}

peptide <- prolfquapp::read_DIANN_output(
  diann.path = diann.path,
  fasta.file = fasta.files,
  nrPeptides = REPORTDATA$nrPeptides,
  Q.Value = 0.01 )

prot_annot <- prolfquapp::dataset_protein_annot(
  peptide,
  c("protein_Id" = "Protein.Group"),
  protein_annot = "fasta.header",
  more_columns = c("nrPeptides", "fasta.id", "Protein.Group.2", "protein_length")
)


annot$raw.file[ !annot$raw.file %in% sort(unique(peptide$raw.file)) ]
nr <- sum(annot$raw.file %in% sort(unique(peptide$raw.file)))
logger::log_info("nr : ", nr, " files annotated")
annot$Relative.Path <- NULL
peptide <- dplyr::inner_join(annot, peptide, multiple = "all")

#### this comes from DIANN
atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("Protein.Group")
atable$hierarchy[["peptide_Id"]] <- c("Stripped.Sequence")
atable$set_response("Peptide.Quantity")
atable$hierarchyDepth <- 1

res <- prolfquapp::dataset_set_factors_deprecated(atable, peptide, SAINT = TRUE)

# attach annotation to combined_protein data
# annotation$raw.file <- basename(annotation$relative.path)
# annotation <- dplyr::mutate(annotation, raw.file = gsub(".raw", "", raw.file))
# annotation$relative.path <- NULL

atable <- res$atable
peptide <- res$msdata

# Preprocess data - aggregate proteins.
config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(peptide, config)

lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()


logger::log_info("AGGREGATING PEPTIDE DATA!")

lfqdataProt <- prolfquapp::aggregate_data(
  lfqdata,
  agg_method = "topN")


lfqdataProt <- prolfquasaint::normalize_exp(lfqdataProt, normalization = REPORTDATA$Normalization)
lfqdataProt <- prolfquasaint::transform_force(lfqdataProt, transformation = REPORTDATA$Transformation)

lfqdata <- lfqdataProt
protAnnot <- prolfqua::ProteinAnnotation$new(lfqdata , prot_annot)

RESULTS <- list() # RESULT is stored in excel table
RESULTS$annotation <- lfqdata$factors()

intdata <- dplyr::inner_join(protAnnot$row_annot, lfqdata$data, multiple = "all")

#debug(protein_2localSaint)
localSAINTinput <- prolfquasaint::protein_2localSaint(
  intdata,
  quantcolumn = lfqdata$config$table$get_response(),
  proteinID = "protein_Id",
  proteinLength = "protein_length",
  IP_name = "raw.file",
  baitCol = "Bait_",
  CorTCol = "CONTROL")

localSAINTinput$inter$exp_transformedIntensity |> summary()

RESULTS <- c(RESULTS, localSAINTinput)
resSaint <- prolfquasaint::runSaint(localSAINTinput, spc = REPORTDATA$spc)


resSaint$list <- dplyr::inner_join(protAnnot$row_annot, resSaint$list,
                                   by = c(protein_Id = "Prey"),
                                   keep = TRUE , multiple = "all")

#resSaint$list$protein_Id <- NULL

RESULTS <- list() # RESULT is stored in excel table
RESULTS$annotation <- lfqdata$factors()

RESULTS <- c(RESULTS, resSaint)
# write analysis results

# Prepare result visualization and render report
cse <- prolfquasaint::ContrastsSAINTexpress$new(resSaint$list)

resContrasts <- cse$get_contrasts()



sig <- resContrasts |>
  dplyr::filter(.data$BFDR  <  REPORTDATA$FDRthreshold &
                  .data$log2_EFCs  >  log2(REPORTDATA$FCthreshold))

# Transform data for PCA visualization etc
tt <- lfqdata$get_Transformer()$log2()
lfqdata_transformed <- tt$lfq

REPORTDATA$pups <- prolfqua::UpSet_interaction_missing_stats(
  lfqdataProt$data,
  lfqdata$config,tr = 2)
RESULTS$InputData <- lfqdata$to_wide()$data

gs <- lfqdata$get_Summariser()
RESULTS$MissingInformation <- gs$interaction_missing_stats()$data
RESULTS$MissingInformation$isotopeLabel <- NULL
RESULTS$listFile <- NULL
writexl::write_xlsx(RESULTS, path = file.path(ZIPDIR,paste0(treat, "_data.xlsx")))

REPORTDATA$BFABRIC <- BFABRIC
REPORTDATA$lfqdata_transformed <- lfqdata_transformed
REPORTDATA$sig <- sig
REPORTDATA$resContrasts <- resContrasts
REPORTDATA$prot_annot <- dplyr::rename(prot_annot, protein = protein_Id)


tmp <- REPORTDATA$pups$data
tmp$UniprotID <- sapply(tmp$protein_Id, function(x){strsplit(x,";")[[1]][1]})

write.table(data.frame(tmp$UniprotID), file = file.path(ZIPDIR,"ORA_background.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE )
sig |> dplyr::group_by(Bait) |> tidyr::nest() -> sigg

if (nrow(sigg) > 0) {
  for (i in 1:nrow(sigg)) {
    tmp <- sigg$data[[i]]
    tmp$UniprotID <- sapply(tmp$Prey, function(x){strsplit(x,";")[[1]][1]})

    filename <- paste0("ORA_Bait_", sigg$Bait[i] , ".txt")
    write.table(data.frame(tmp$UniprotID),
                file = file.path(ZIPDIR, filename),
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE )
  }
}

prolfquasaint::copy_SAINTe_doc(workdir = ZIPDIR)

SEP <- REPORTDATA

saveRDS(REPORTDATA, file = "REPORTDATA.rds")
rm(list = setdiff(ls(), c("REPORTDATA","ZIPDIR","treat","yml"))) # delete all variables not needed for rendering
SEP <- REPORTDATA

fragPipeDIA <- names(yml$application$input) == "FragPipeDIA-dataset"
text <-  c( "The LC-MS data was processed using the ",
            if(fragPipeDIA){"FragPipe-DIA proteomics pipeline [@yu2023analysis]."} else {"DIA-NN software [@demichev2020dia]."} ,
            "The quantification results were extracted from the DIANN main report, containing precursor (charged modified peptide) quantification results.",
            "The protein abundances were estimated by the sum of all the precursor abundances (Precursor.Quantity column) assigned to a protein. ")

text <- paste(text, collapse = "")

rmarkdown::render("SaintExpressReportMsFragger.Rmd",
                  params = list(sep = REPORTDATA, textpreprocessing = text),
                  output_format = bookdown::html_document2())


file.copy("SaintExpressReportMsFragger.html",
          file.path(ZIPDIR, paste0(treat, "SaintExpressReportMsFragger.html")),
          overwrite = TRUE)

