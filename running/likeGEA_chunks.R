library(likeGEA)

base <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/angsd/snp_calling_global"

beagle_file <- file.path(
  base,
  "current.combined.subsetted.norep.unrel.noPP.noLP.beagle.gz"
)

chunk_dir <- file.path(base, "likeGEA_chunks_4")
temp_dir <- file.path(chunk_dir, "temp_chunks")

stopifnot(file.exists(beagle_file))

dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

# Start with a clean temporary directory
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
}

dir.create(temp_dir)

old_wd <- getwd()
setwd(temp_dir)

chunks <- get_chunks(
  beagle_file = beagle_file,
  loci_per_chunk = 25000
)

setwd(old_wd)

# Get generated files and sort numerically
chunks <- list.files(
  temp_dir,
  pattern = "^chunk_[0-9]+\\.beagle\\.gz$",
  full.names = TRUE
)

chunk_number <- as.integer(
  sub("^chunk_([0-9]+)\\.beagle\\.gz$", "\\1", basename(chunks))
)

chunks <- chunks[order(chunk_number)]

stopifnot(length(chunks) > 0)

# -----------------------------------
# Standardize chunk filenames

# Only replace existing chunks after the
# new chunking operation has succeeded.
unlink(
  list.files(
    chunk_dir,
    pattern = "^chunk_[0-9]+\\.beagle\\.gz$",
    full.names = TRUE
  )
)

for (i in seq_along(chunks)) {
  
  new_file <- file.path(
    chunk_dir,
    sprintf("chunk_%03d.beagle.gz", i - 1)
  )
  
  if (!file.rename(chunks[i], new_file)) {
    stop("Could not move chunk to: ", new_file)
  }
}

unlink(temp_dir, recursive = TRUE)

message("Created ", length(chunks), " chunks.")
message("Saved in: ", chunk_dir)