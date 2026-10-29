#' Format a beagle file as phred-scaled likelihoods.
#'
#' Takes a beagle formatted input file and converts it to a matrix of phred-
#' scaled likelihoods, as needed by \code{mechGEA}. See details.
#'
#' The input file format required by this package has one row for each locus and
#' one column for each individual. Phred-scaled likelihoods (PL) for the three
#' possible genotypes at each locus/individual are comma seperated (major
#' homozygote, heterozygote, minor heterozygote). A higher PL constitutes a
#' lower likelihood; by convention the minimum PL is subtracted from each
#' triplet such that the most likely genotype has a PL of 0.
#'
#' Phred-likelihoods can be calculated according to:
#' \deqn{PL(\text{genotype}) = -10 \times \log_{10}\left( P(\text{Data} \mid \text{Genotype}) \right)}
#'
#' Where \eqn{P(\text{Data} \mid \text{Genotype})} is the a genotype likelihood
#' as produced by something like \code{ANGSD}'s \code{doGlf}.
#'
#' @param file character. File path to the beagle or beagle.gz file containing
#'   the raw data.
#' @param lines numeric, default NULL. Optional vector of numbers noting which
#'   lines in the beagle file to keep. This must be sequential due to
#'   \code{\link[data.table]{fread}}'s restrictions.
#' @param ind_order numeric, default NULL. Numeric vector indicating final
#'   sorted order to use for individuals the output matrix (\code{c(2,1,3)}, for
#'   example, would put the second individual in \code{file} first, then the
#'   first individual second, and then the third individual last in the final
#'   matrix, for example). If NULL, the order will be unchanged.
#' @param ind_names character, default NULL. Names for individuals in the final
#'   matrix. If \code{ind_order} is provided, this should match the initial
#'   order of sample names, not the final order.
#'
#' @return A \eqn{p \times n} matrix, where p is the number of loci and n is
#'   the number of individuals. See details for format.
#'
#' @export
#'
#' @author William Hemstrom
#'
#' @examples
#' # example code
#'
#' beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
#' pl <- function(beagle,
#'                lines = 1:100) # first 100 lines
#'
beagle2pl <- function(file, lines = NULL, ind_order = NULL, ind_names = NULL){

  if(!is.null(lines)){
    d <- fread(file, skip = lines[1] - 1, nrows = length(lines))
  }
  else{
    d <- fread(file)
  }

  header <- d[,1:3]
  d <- d[,-c(1:3)]
  d <- as.matrix(d)
  d[d == 0] <- 0.00000000001
  d <- -10*log10(d)
  d <- round(d)
  storage.mode(d) <- "integer"

  # loop through inds to standardize, should be quickish
  inds <- 1:ncol(d)
  inds <- split(inds, rep(1:(ncol(d)/3), each = 3))
  for(i in 1:length(inds)){
    m <- matrixStats::rowMins(d[,inds[[i]]])
    d[,inds[[i]]] <- d[,inds[[i]]] - m
  }

  # convert to proper strings
  storage.mode(d) <- "character"
  dfinal <- matrix(as.character(NA), nrow(d), ncol(d)/3)
  for(i in 1:length(inds)){
    dfinal[,i] <- paste0(d[,inds[[i]][1]], ",", d[,inds[[i]][2]], ",", d[,inds[[i]][3]])
  }

  dfinal[dfinal == "0,0,0"] <- "."

  # names
  if(!is.null(ind_names)){
    colnames(dfinal) <- ind_names
  }

  # reorganize if requested
  if(!is.null(ind_order)){
    dfinal <- dfinal[,ind_order]
  }

  return(dfinal)
}

#' Chunk a beagle file into partitions for processing.
#'
#' Creates n chunks of a beagle file, where sequential lines are sent to
#' different files round-robin style. Note that this is much faster on Linux
#' than on Windows, since it avoids reading anything into memory using unix
#' commands.
#'
#' @param beagle_file path to \code{.beagle} or \code{.beagle.gz} file.
#' @param loci_per_chunk numeric, default 10000. Number of loci to distribute to each
#'   chunk.
#' @param max_read numeric, default 50000. Number of lines to read for splitting
#'   at once. Ignored on Linux.
#'
#' @return A vector of filenames for the chunks.
#'
#' @export
#'
#' @author William Hemstrom
get_chunks <- function(beagle_file,
                       loci_per_chunk = 10000,
                       max_read = 50000){

  sys <- Sys.info()["sysname"]

  # get number of lines
  if(tools::file_ext(beagle_file) == "gz"){
    if(sys == "Linux"){
      nsnps <- system(paste0("zcat ", beagle_file, " | wc -l"), intern = TRUE)
      nsnps <- as.numeric(nsnps) - 1
    }
    else {
      nsnps <- fread(beagle_file, select = 1)
      nsnps <- nrow(nsnps)
    }
  } else {
    nsnps <- bread::bnrow(beagle_file)
    nsnps <- nsnps - 1
  }

  if(nsnps <= loci_per_chunk){
    cat("File has less loci than requested per-chunk already. Nothing to do.\n")
    return(beagle_file)
  }

  #===========split==========

  nfiles <- ceiling(nsnps/loci_per_chunk)

  cat("Splitting into", nfiles, "files with approx.", nsnps/nfiles, "snps per file.\n")

  # Linux -- use cmd
  if(sys == "Linux"){
    if(tools::file_ext(beagle_file) == "gz"){
      cmd <- paste0(
        "f=", shQuote(beagle_file), "; ",
        "export HEADER=$(zcat \"$f\" 2>/dev/null | head -n 1); ",
        "zcat \"$f\" | tail -n +2 | ",
        "split -n r/", nfiles, " -d --additional-suffix=.beagle.gz ",
        "--filter='{ printf \"%s\\n\" \"$HEADER\"; cat; } | gzip > $FILE' ",
        "- chunk_"
      )
    } else {
      cmd <- paste0("tail -n +2 ", shQuote(beagle_file),
                    " | split -n r/", nfiles,
                    " -d --additional-suffix=.beagle --filter='{ head -n 1 ",
                    shQuote(beagle_file), " ; cat; } > $FILE' - chunk_")
    }

    system(cmd)

    outfiles <- paste0("chunk_", formatC(0:(nfiles - 1), width = nchar(as.character(nfiles)), format = "d", flag = "0"),
                       ".beagle")
    if(tools::file_ext(beagle_file) == "gz"){outfiles <- paste0(outfiles, ".gz")}
    return(outfiles)
  }

  # windows is sucky, read in chunks and re-write since we don't have access to split.
  # NOTE if Nikunj confirms that there isn't a downside to doing this in sequential chunks,
  # then we can unzip, then use bread to split very quickly and then rezip.
  else if(sys == "Windows"){
    fnames <- paste0("chunk_", 1:nfiles, ".beagle.gz")

    read_parts <- split(1:nsnps, sort(rep(1:(nsnps/max_read), length.out = nsnps)))
    read_parts <- lapply(read_parts, function(x){
      p <- c(min(x), max(x))
      p <- c(p[1], p[2] - p[1])
      return(p)
    })
    for(i in 1:length(read_parts)){
      cat("Reading part:", i)
      tbeag <- fread(beagle_file, skip = read_parts[[i]][1] - 1,
                     nrows = read_parts[[i]][2])

      cparts <- split(tbeag, rep(1:nfiles, length.out = nrow(tbeag)))
      names(cparts) <- fnames

      lapply(1:length(cparts), function(x){

        fwrite(cparts[[x]], names(cparts)[x],
               append = ifelse(file.exists(names(cparts)[x]), TRUE, FALSE),
               sep = "\t")
      })
    }
    return(fnames)
  }
}

#' Format metadata for analysis.
#'
#' From a \code{\link[sf]{sf}} object containing patch and individual IDs as
#' well as environmental variables, create the three metadata input files
#' required by \code{\link{mechGEA}}. An example input is in
#' \code{\link{sample_metadata}}. An example output is in
#' \code{\link{example_gl}}.
#'
#' @param sample_meta \code{\link[sf]{sf}} POINT geometry containing two
#'   required columns: \code{patch} containing patch IDs and \code{ind_id}
#'   containing individual IDs. Additional columns specify named environmental
#'   variables. Patches must be in continuous blocks and sorted identically to
#'   the likelihood data also provided to \code{\link{mechGEA}}. NOTE:
#'   \code{.beagle} or \code{.beagle.gz} files can be re-sorted on import using
#'   the \code{\link{beagle2pl}} function.
#' @param env_mean_method character, default "unweighted". Specifies how
#'   environmental variables are averaged across samples. Options:
#'   * "unweighted": direct means across samples within patches.
#'   * "weighted": inverse-distance scaled averages. Scaling power given by
#'     \code{idp}.
#' @param idp integer, default 1. Inverse-distance power used for weighted
#'   environmental averages. Ignored if \code{env_mean_method == "weighted"}.
#'   Higher values give less weight to points distant from the patch centroid.
#'
#' @export
#' @author William Hemstrom
prepare_metadata <- function(sample_meta,
                             env_mean_method = "unweighted",
                             idp = 1){

  patch <- ind_id <- lat <- lon <- ..targets <- NULL

  # first things first--check the sorting and if it isn't OK. We don't want to re-sort,
  # since the genotypes should be in the same order and that could cause issues

  # patch info per sample
  dt <- as.data.table(sample_meta)
  patch_info <- dt[,.(ind_id, patch)]
  in_blocks <- TRUE
  observed <- character()
  for(i in 1:nrow(patch_info)){
    tpatch <- patch_info$patch[i]
    if(tpatch %notin% observed){
      observed <- c(observed, tpatch)
      next
    } else {
      if(any(observed[-length(observed)] == tpatch)){
        in_blocks <- FALSE
      }
    }
  }

  if(isFALSE(in_blocks)){
    stop("Some observations are not in continous blocks by patch. Individuals in both the metadata AND the genotypic
         data MUST be in blocks/sorted by patch.\n")
  }



  # table of patches with x/y, sorted correctly
  # Coordinates
  patch_xy <- data.table(patch = unique(sample_meta$patch),
                         lat = 0,
                         lon = 0)
  for(i in 1:nrow(patch_xy)){
    tm <- sf::st_centroid(sf::st_union(sample_meta[sample_meta$patch == patch_xy$patch[i],]))
    tm <- sf::st_coordinates(tm)
    patch_xy[i,lat := tm[1,2]]
    patch_xy[i,lon := tm[1,1]]
  }

  # mean envs per patch
  targets <- colnames(dt)
  targets <- targets[targets %notin% c("patch", "ind_id", "geometry")]
  if(env_mean_method == "unweighted"){
    envs <- dt[,lapply(.SD, mean), by = patch, .SDcols = targets]
    ptchs <- envs$patch
    targets <- c("patch", targets)
    .fix..call(envs <- envs[,..targets])
  } else if(env_mean_method == "weighted"){ # for weighted (idw); relies on gstat
    envs <- copy(patch_xy)[,1]
    for(i in targets) envs[[i]] <- 0

    # loop through and get dist to centroid for each point
    for(i in 1:nrow(patch_xy)){
      this_sf <- sample_meta[sample_meta$patch == patch_xy$patch[i],]
      this_center <- sf::st_as_sf(patch_xy[i,], coords = c("lon", "lat"), crs = sf::st_crs(sample_meta))

      for(j in 1:length(targets)){
        idw_result <- gstat::idw(
          formula = as.formula(paste0(targets[j], " ~ 1")),
          locations = this_sf,
          newdata = this_center,
          idp = idp)[[1]]

        set(envs, i = i, j = targets[[j]], value = idw_result)
      }
    }
  } else {
    stop("Unrecognized env_mean_method. Options: 'weighted', 'unweighted'.\n")
  }

  # return
  ret <- list(patch_info = patch_info[,.(ind_id, patch)],
              patches_xy = patch_xy,
              env = envs)
  return(lapply(ret, as.data.frame))
}


