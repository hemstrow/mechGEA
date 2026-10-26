test_that("get_chunks",{

  # correct on linux with gz
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  chunk_files <- get_chunks(beagle, 100)

  ## ref
  ref <- fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, fread)
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
  ref <- fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, fread)
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
