# from snpR
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

# from snpR
.fix..call <- function(fun) {
  return(.suppress_specific_warning(fun, "variable in calling scope for clarity"))
}

#' Check that the mechGEA package is fully, correctly installed.
#'
#' This function runs a series of checks to see that mechGEA is correctly
#' installed. It runs automatically during the \code{\link{mechGEA}} pipeline,
#' but can also be run manually. It only needs to be run once per session.
#'
#' @param status install status, fetched from a system environment (which should
#'   not be manually updated; this function will check and update the status as
#'   needed).
#'
#' @author William Hemstrom
#' @export
check_correct_install <- function(status = Sys.getenv("mechGEA_full_install")){
  if(status == ""){
    # do a quick, tiny test compilation under rstan and stan to figure out what runs. Most functions are turned-off unless these are enabled.
    ## rstan; can we expose a function?
    Sys.setenv(mechGEA_full_install = TRUE)

    if(!"rstan" %in% rownames(utils::installed.packages())){
      Sys.setenv(mechGEA_full_install = FALSE)
      warning("rstan not installed. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
    } else {
      .calc_sel_code <- "
    functions{
      real x_plus_y(real x, real y){
        real xysum;
        xysum = x + y;
        return xysum;
      }
    }
    "
      OK <- try(rstan::expose_stan_functions(rstan::stanc(model_code = .calc_sel_code)), silent = TRUE)
      if(methods::is(OK, "try-error") | !exists("x_plus_y")){
        Sys.setenv(mechGEA_full_install = FALSE)
        warning("rstan installed but not working properly. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
      }
    }


    ## cmdstanr
    if(!"cmdstanr" %in% rownames(utils::installed.packages())){
      Sys.setenv(mechGEA_full_install = FALSE)
      warning("cmdstanr not installed. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
    } else {
      OK <- try(cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE), silent = TRUE)
      if(!is.null(OK) | methods::is(OK, "try-error")){
        Sys.setenv(mechGEA_full_install = FALSE)
        warning("cmdstanr could not find cmdstan toolchain. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
      } else {
        OK <- try(cmdstanr::cmdstan_path(), silent = TRUE)
        if(methods::is(OK, "try-error")){
          Sys.setenv(mechGEA_full_install = FALSE)
          warning("cmdstanr could not locate cmdstan. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
        } else {
          file <- file.path(cmdstanr::cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
          mod <- try(cmdstanr::cmdstan_model(file), silent = TRUE)
          if(methods::is(mod, "try-error")){
            Sys.setenv(mechGEA_full_install = FALSE)
            warning("cmdstanr could not build a test model. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
          } else {
            data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))
            fit <- try(
              mod$sample(
                data = data_list,
                seed = 123,
                chains = 2,
                parallel_chains = 2,
                show_messages = FALSE,
                refresh = 0),
              silent = FALSE)
            if(methods::is(fit, "try-error")){
              Sys.setenv(mechGEA_full_install = FALSE)
              warning("cmdstanr could not fit a test model. Many features will be disabled.\nConsider using setup_mechGEA_env().\n")
            }
          }
        }
      }
    }

    status <- Sys.getenv("mechGEA_full_install")
  }

  if(as.logical(status)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# from snpR
.check.installed <- function(pkg, install.type = "basic", source = NULL){
  for(i in 1:length(pkg)){

    tpkg <- pkg[i]

    if(!tpkg %in% rownames(utils::installed.packages())){
      say <- paste0("Package '", tpkg, "' not found.")
      cat(say, "")
      cat("Install? (y or n)\n")

      if(!interactive()){
        stop(say)
      }

      resp <- readLines(n = 1)
      resp <- tolower(resp)


      if(resp != "y"){
        stop(say)
      }
      else{
        utils::install.packages(tpkg)
      }
    }
  }

  return(TRUE)
}

