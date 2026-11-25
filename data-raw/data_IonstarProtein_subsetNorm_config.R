# Rebuild the configuration of data_IonstarProtein_subsetNorm with the current AnalysisConfiguration class
#
# Problem: the stored configuration was serialized by an old prolfqua version (camelCase fields and
# `ident_q_value = "pep"`, a column the protein data does not have). prolfqua 1.9.11 replaced `ident_q_value` by
# `identification_confidence`.
#
# Fix: build a fresh configuration from the old field values. The protein abundances carry no identification
# evidence, so no identification confidence column is added.
# Run from the package root.

load("data/data_IonstarProtein_subsetNorm.rda")
old <- data_IonstarProtein_subsetNorm$config

config <- prolfqua::AnalysisConfiguration$new()
config$file_name <- old$file_name
config$sample_name <- old$sample_name
config$hierarchy <- old$hierarchy
config$hierarchy_depth <- old$hierarchy_depth
config$factors <- old$factors
config$factor_depth <- old$factor_depth
config$work_intensity <- old$work_intensity
config$is_response_transformed <- old$is_response_transformed
config$nr_children <- old$nr_children
config$sep <- old$sep

data_IonstarProtein_subsetNorm$config <- config
save(data_IonstarProtein_subsetNorm, file = "data/data_IonstarProtein_subsetNorm.rda", compress = "xz")
