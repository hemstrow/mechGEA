test_that("get_chunks",{

  # correct on linux with gz
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  chunk_files <- get_chunks(beagle, 100)

  ## ref
  ref <- data.table::fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, data.table::fread)
  ## colnames
  for(x in 1:length(chunks)){
    expect_identical(colnames(chunks[[x]]), colnames(ref))
  }
  ## contents
  interleaved <- rbindlist(chunks, idcol = "chunk")[, row_num := seq_len(.N), by = chunk]
  interleaved <- interleaved[order(row_num, chunk)]
  interleaved$chunk <- interleaved$row_num <- NULL
  expect_identical(interleaved, ref)
  file.remove(chunk_files)

  # correct on linux with non-gz
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  file.copy(beagle, "test.beagle.gz")
  system("gunzip test.beagle.gz")
  beagle <- "test.beagle"
  chunk_files <- get_chunks(beagle, 100)

  ## ref
  ref <- data.table::fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, data.table::fread)
  ## colnames
  for(x in 1:length(chunks)){
    expect_identical(colnames(chunks[[x]]), colnames(ref))
  }
  ## contents
  interleaved <- rbindlist(chunks, idcol = "chunk")[, row_num := seq_len(.N), by = chunk]
  interleaved <- interleaved[order(row_num, chunk)]
  interleaved$chunk <- interleaved$row_num <- NULL
  expect_identical(interleaved, ref)
  file.remove(chunk_files)
  file.remove("test.beagle")
  if(file.exists("test.beagle.gz")){file.remove("test.beagle.gz")}
})

test_that("env pruning",{
  set.seed(123122)
  # generate correlated envs
  X <- rnorm(100)

  # Create Y as a linear combination to inject correlation
  r <- rnorm(40, ifelse(rbinom(20, 1, .5) == 1, .9, -.7), .2)
  r[r >= 1] <- .99999
  r[r <= -1] <- -.99999
  Y <- matrix(0, 100, length(r) + 1)
  Y[,1] <- X
  for(i in 2:ncol(Y)){
    Y[,i] <- r[i-1] * Y[,i-1] + sqrt(1 - r[i-1]^2) * rnorm(100)
  }
  Y <- as.data.table(Y)
  colnames(Y) <- paste0("env", 1:ncol(Y))


  # test
  ty <- cor(Y)
  ty[ty == 1] <- 0
  expect_true(max(abs(ty)) >= .6)
  res <- prune_envs(Y, r2 = .7, seed = 42)
  tx <- cor(res)
  tx[tx == 1] <- 0
  expect_true(max(abs(tx)) <= .7)

  # test with a patch column -- same result given same seed, but with a patch column
  Ypatch <- cbind(patch = paste0("patch", 1:nrow(Y)), Y)
  res2 <- prune_envs(Ypatch, r2 = .7, seed = 42)
  expect_identical(res, res2[,-1])
  expect_true(colnames(res2)[1] == "patch")


})

test_that("beagle2pl",{

  # since the example beagle was made by converting the pl example into a beagle, we can do a return-trip test for formatting.
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  pl <- beagle2pl(beagle, ind_names = colnames(example_gl$likelihoods))

  expect_equal(pl, example_gl$likelihoods, ignore_attr = TRUE)
  expect_equal(colnames(pl), colnames(example_gl$likelihoods))

  # test nline selection
  pl <- beagle2pl(beagle, ind_names = colnames(example_gl$likelihoods), lines = 1:100)
  expect_equal(pl, example_gl$likelihoods[1:100,], ignore_attr = TRUE)
  expect_equal(colnames(pl), colnames(example_gl$likelihoods))

  # test re-ordering
  pl <- beagle2pl(beagle, ind_names = colnames(example_gl$likelihoods), ind_order = c(101:200, 1:100))
  expect_equal(pl, example_gl$likelihoods[,c(101:200, 1:100)], ignore_attr = TRUE)
  expect_equal(colnames(pl), colnames(example_gl$likelihoods)[c(101:200, 1:100)])

})

test_that("prepare_metadata",{

  # test that the basic formatting is OK
  input <- prepare_metadata(sample_metadata)

  expect_equal(input$patch_info, example_gl$patch_info)
  expect_equal(input$patches_xy, example_gl$patches_xy)
  expect_equal(input$env, example_gl$env)

  # test that unsorted data throws an error
  test <- sample_metadata[sample(nrow(sample_metadata), nrow(sample_metadata), FALSE),]
  expect_error(prepare_metadata(test), "continous blocks by patch")

  # weighted vs unweighted
  ## generate varaince in coords
  mx <- c(-160, -165)
  my <- c(45, 50)
  n <- nrow(sample_metadata)
  set.seed(42)
  ndt <- data.table(patch = c(rep("patch_001", n/2),
                              rep("patch_002", n/2)),
                    ind_id = sample_metadata$ind_id,
                    lon = c(rnorm(n/2, mx[1], .5),
                            rnorm(n/2, mx[2], .5)),
                    lat = c(rnorm(n/2, my[1], .5),
                            rnorm(n/2, my[2], .5)),
                    env1 = rnorm(n, 0, 1),
                    env2 = rnorm(n, 0, 1))
  ndt <- sf::st_as_sf(ndt, coords = c("lon", "lat"), crs = sf::st_crs(sample_metadata))

  uw <- prepare_metadata(ndt)
  expect_equal(uw$env$env1, c(mean(ndt[ndt$patch == "patch_001",]$env1),
                              mean(ndt[ndt$patch == "patch_002",]$env1)))
  expect_equal(uw$env$env2, c(mean(ndt[ndt$patch == "patch_001",]$env2),
                              mean(ndt[ndt$patch == "patch_002",]$env2)))
  w1 <- prepare_metadata(ndt, env_mean_method = "weighted", idp = 1)
  w4 <- prepare_metadata(ndt, env_mean_method = "weighted", idp = 4)
  w0 <- prepare_metadata(ndt, env_mean_method = "weighted", idp = 0)
  d1 <- sum((uw$env[,-1] - w1$env[,-1])^2)
  d4 <- sum((uw$env[,-1] - w4$env[,-1])^2)
  expect_true(d1 < d4) # lower idp should always be closer to straight means
  expect_equal(uw$env, w0$env) # with zero weight, will be the mean

})
