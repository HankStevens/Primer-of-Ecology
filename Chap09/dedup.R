### deduplicate citations
library(synthesisr)

imported_files <- read_refs(
  filename = "~/trial.bib",
  return_df = TRUE)

system.time(
  dups <- find_duplicates(
    imported_files$title,
    method = "string_osa",
    rm_punctuation = TRUE,
    to_lower = TRUE
  )
)
system.time(
  clean <- deduplicate(
    imported_files,
    match_by="title",
    method = "string_osa",
    rm_punctuation = TRUE,
    to_lower = TRUE
  )
)
write_refs(clean, format="bib", file="trialClean.bib")
