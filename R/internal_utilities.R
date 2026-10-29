# checks if cmdstanr can find cmdstan and stops if not. Otherwise returns path.
.check_cmdstan <- function(){
  path <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  if(methods::is(path, "try-error")){

    say <- "No cmdstan install found.\n"

    if(!interactive()){
      stop(say)
    }

    cat(say, "")
    cat("Install? (y/n)\n")

    resp <- readLines(n = 1)
    resp <- tolower(resp)
    if(resp == "yes") resp <- "y"
    if(resp != "y"){
      stop(say)
    } else{
      cmdstanr::install_cmdstan()
    }

    path <- try(cmdstanr::cmdstan_version(), silent = TRUE)
    if(methods::is(path, "try-error")){
      stop("cmdstan still cannot be located. Likely installation failure. Try a manual install?")
    }
  }
  return(path)
}

.console_hline <- function(char = "="){
  return(paste0(rep(char, getOption("width")), collapse = ""))
}


# from snpR
.suppress_specific_warning <- function(fun, warn_to_suppress) {
  warn_handle <- function(warn) {
    if (any(grepl(warn_to_suppress, warn))) {
      invokeRestart("muffleWarning")
    }
  }
  withCallingHandlers(res <- fun, warning = warn_handle)
  return(res)
}
.fix..call <- function(fun) {
  return(.suppress_specific_warning(fun, "variable in calling scope for clarity"))
}
