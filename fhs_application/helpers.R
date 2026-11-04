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
