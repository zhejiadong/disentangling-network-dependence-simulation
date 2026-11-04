# Offspring Exam 2: raw phenotype, common complete cohort, then two networks.
# Work from fhs_application/. All generated objects stay in the current environment.
library(igraph)
# Selective dbGaP table reads and person-ID validation.
read_dbgap_columns <- function(path, columns) {
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  text <- readLines(con, n = 60L, warn = FALSE)
  close(con)
  header <- which(vapply(strsplit(text, "\t"), function(x) "shareid" %in% x, logical(1)))[1L]
  if (is.na(header)) stop("No shareid header in ",basename(path))
  all_names <- strsplit(text[header], "\t")[[1L]]
  if (!all(columns %in% all_names)) stop("Missing required columns in ",basename(path))
  classes <- rep("NULL", length(all_names))
  classes[match(columns, all_names)] <- "character"
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  read.delim(con, skip = header - 1L, colClasses = classes, quote = "",
             comment.char = "", check.names = FALSE, na.strings = c("", ".", "NA"),
             stringsAsFactors = FALSE)
}

assert_unique_shareids <- function(data, label) {
  data$shareid <- trimws(as.character(data$shareid))
  if (anyNA(data$shareid) || any(!nzchar(data$shareid))) stop(label, " contains missing IDs.")
  if (anyDuplicated(data$shareid)) stop(label, " contains duplicate IDs after c1/c2 bind.")
  data
}

read_consent_columns <- function(paths, columns) {
  pieces <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    pieces[[i]] <- read_dbgap_columns(paths[i], columns)
    pieces[[i]]$consent_group <- names(paths)[i]
  }
  do.call(rbind, pieces)
}

# Put the authorized c1/c2 dbGaP tables (plain text or .gz) under data_raw/.
raw_paths <- list.files("data_raw", recursive = TRUE, full.names = TRUE)
table_paths <- raw_paths[!grepl("[.]xml$", raw_paths, ignore.case = TRUE)]
files <- list()
for (dataset in c("pht000031", "pht003099", "pht000836")) {
  dataset_paths <- table_paths[grepl(dataset, basename(table_paths))]
  files[[dataset]] <- c(
    c1 = dataset_paths[grepl("(^|[._-])c1([._-]|$)", basename(dataset_paths))],
    c2 = dataset_paths[grepl("(^|[._-])c2([._-]|$)", basename(dataset_paths))]
  )
  if (length(files[[dataset]]) != 2L || !identical(names(files[[dataset]]),c("c1","c2"))) {
    stop("Supply exactly one c1 and one c2 table for ",dataset)
  }
}
exam2 <- read_consent_columns(files$pht000031,
  c("shareid", "B13", "B14", "B24", "B439", "B443", "B448"))
demographics <- read_consent_columns(files$pht003099, c("shareid", "age2", "sex"))

exam2 <- assert_unique_shareids(exam2, "Exam-2 phenotype")
demographics <- assert_unique_shareids(demographics, "Demographics")

# BMI in kg/m^2; LVM/height in g/m; SBP is first-physician B24.
phenotype <- data.frame(
  id = trimws(as.character(exam2$shareid)),
  consent_exam2 = exam2$consent_group,
  weight_kg = as.numeric(exam2$B13),
  height = as.numeric(exam2$B14) / 100,
  sbp = as.numeric(exam2$B24),
  septal_thickness_mm = as.numeric(exam2$B439),
  posterior_wall_thickness_mm = as.numeric(exam2$B443),
  lvid_mm = as.numeric(exam2$B448),
  stringsAsFactors = FALSE
)
phenotype$bmi <- phenotype$weight_kg / phenotype$height^2
phenotype$lvm_raw <- 1.04 * (
  (phenotype$lvid_mm + phenotype$septal_thickness_mm +
     phenotype$posterior_wall_thickness_mm)^3 - phenotype$lvid_mm^3
) / 1000 - 13.6
phenotype$lvm_height1 <- phenotype$lvm_raw / phenotype$height
for (variable in c("height", "bmi", "lvm_raw", "lvm_height1", "sbp")) {
  phenotype[[variable]][!is.finite(phenotype[[variable]])] <- NA_real_
}
demo <- data.frame(
  id = trimws(as.character(demographics$shareid)),
  consent_demo = demographics$consent_group,
  age = as.numeric(demographics$age2),
  sex = as.numeric(demographics$sex),
  stringsAsFactors = FALSE
)
joined <- merge(phenotype, demo, by = "id", all = FALSE, sort = FALSE)

stopifnot(!anyDuplicated(joined$id), all(joined$consent_exam2 == joined$consent_demo))

# Height is required to construct BMI and the indexed outcome, not an extra covariate.
required <- c("id", "bmi", "lvm_raw", "lvm_height1", "age", "sex", "height", "sbp")
complete <- complete.cases(joined[, required])
cohort_all <- joined[complete, c("id", "consent_exam2", "consent_demo", "bmi",
  "lvm_raw", "lvm_height1", "age", "sex", "height", "sbp")]
rownames(cohort_all) <- NULL
if (!nrow(cohort_all)) stop("The maximal complete P1 cohort is empty.")
if (any(!cohort_all$sex %in% c(1, 2))) stop("Sex contains values outside 1/2.")
if (any(cohort_all$height < 0.9 | cohort_all$height > 2.3)) stop("Height outside 0.9-2.3 m.")
if (any(cohort_all$bmi < 10 | cohort_all$bmi > 80)) stop("BMI outside 10-80 kg/m^2.")
if (any(cohort_all$lvm_raw <= 0) || any(cohort_all$lvm_height1 <= 0)) stop("Invalid LVM outcome.")
if (any(cohort_all$sbp < 60 | cohort_all$sbp > 260)) stop("B24 SBP outside 60-260 mm Hg.")
cohort_all$female <- as.integer(cohort_all$sex == 2L)
cohort_flow <- data.frame(
  stage = c("Exam-2 phenotype rows", "Demographic rows", "Joined rows", "Maximal complete rows"),
  n = c(nrow(phenotype), nrow(demo), nrow(joined), nrow(cohort_all)))
print(cohort_flow, row.names = FALSE)

# Repeated egos are intentional: this is a tie-spell table, not a person table.
raw_network <- read_consent_columns(files$pht000836,
  c("shareid", "sharealterid", "ALTERTYPE", "SPELLBEGIN", "SPELLEND"))
network_rows <- data.frame(
  ego = trimws(as.character(raw_network$shareid)),
  alter = trimws(as.character(raw_network$sharealterid)),
  relationship_type = trimws(as.character(raw_network$ALTERTYPE)),
  spell_begin = as.numeric(raw_network$SPELLBEGIN),
  spell_end = as.numeric(raw_network$SPELLEND))
relationship_types <- list(social = c("FRIENDNR", "COWORKERNR"),
  relatives = c("MOTHER", "FATHER", "CHILD", "BROTHER", "SISTER"))
eligible_rows <- edges <- cohorts <- graphs <- networks <- adjacency <- list()
for (label in names(relationship_types)) {
  eligible <- network_rows$relationship_type %in% relationship_types[[label]] &
    is.finite(network_rows$spell_begin) & is.finite(network_rows$spell_end) &
    network_rows$spell_begin <= 144 & network_rows$spell_end > 108 &
    network_rows$ego %in% cohort_all$id & network_rows$alter %in% cohort_all$id &
    nzchar(network_rows$ego) & nzchar(network_rows$alter) &
    network_rows$ego != network_rows$alter
  if (!any(eligible)) stop("The ", label, " network has zero eligible raw rows.")
  eligible_rows[[label]] <- network_rows[eligible, ]
  edges[[label]] <- unique(data.frame(
    from = pmin(network_rows$ego[eligible], network_rows$alter[eligible]),
    to = pmax(network_rows$ego[eligible], network_rows$alter[eligible])))
  edges[[label]] <- edges[[label]][order(edges[[label]]$from, edges[[label]]$to), ]
  endpoint_ids <- unique(c(edges[[label]]$from, edges[[label]]$to))
  ids <- cohort_all$id[cohort_all$id %in% endpoint_ids]
  cohorts[[label]] <- cohort_all[match(ids, cohort_all$id), ]
  graphs[[label]] <- simplify(graph_from_data_frame(edges[[label]], directed = FALSE,
    vertices = data.frame(name = ids)), remove.multiple = TRUE, remove.loops = TRUE)
  adjacency[[label]] <- as.matrix(as_adjacency_matrix(graphs[[label]], sparse = FALSE, attr = NULL))
  networks[[label]] <- list(ids = ids, graph = graphs[[label]],
    relationship_types = relationship_types[[label]], undirected_edges = edges[[label]])
}
social_cohort <- cohorts$social
relatives_cohort <- cohorts$relatives
social_graph <- graphs$social
relatives_graph <- graphs$relatives
social_adjacency <- adjacency$social
relatives_adjacency <- adjacency$relatives

for (label in names(cohorts)) {
  dir.create(file.path("data", label), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cohorts[[label]], file.path("data", label, "cohort_private.rds"))
  saveRDS(networks[[label]], file.path("data", label, "network_private.rds"))
}
write.csv(cohort_flow, "data/cohort_flow_private.csv", row.names = FALSE)
