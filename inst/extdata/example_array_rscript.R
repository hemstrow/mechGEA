library(mechGEA); library(data.table); library(sf)
args <- commandArgs(TRUE)
beagle <- as.character(args[1])
meta <- readRDS(as.character(args[2]))
parms <- readRDS(as.character(args[3]))
migration_res <- as.character(args[4])

cat("Formatting", beagle, "as phred-scaled likelihoods.\n")

Phred <- beagle2pl(beagle)

cat("Done. Data inputs:\n")
cat("=========================================\n.Phred:\n")
print(utils::head(Phred)[,1:6])
cat("=========================================\n.patch_info:\n")
print(meta$patch_info)
cat("=========================================\n.patches_xy:\n")
print(meta$patches_xy)
cat("=========================================\n.env:\n")
print(meta$env)

use_migration <- FALSE
if(file.exists(migration_res)){
  cat("=========================================\n.migration_res:\n")
  migration_res <- readRDS(migration_res)
  print(migration_res)
  use_migration <- TRUE
}

cat("\n============================\n")
if(length(parms) == 0){
  cat("No additional arguments. Beginning run.\n")
} else {
  cat("Additional arguments:"){
    for(i in 1:length(parms)){
      cat(paste0("\n\t", names(parms)[i], " = ", parms[[i]]))
    }
  }
  cat("\nBeginning run.\n")
}


like_args <- c(list(Phred = Phred,
                    patch_info = meta$patch_info,
                    patches_xy = meta$patches_xy,
                    env = meta$env),
               parms)
if(use_migration){
  like_args <- c(like_args, list(migration_results = migration_res))
}

results <- do.call(mechGEA, like_args)
