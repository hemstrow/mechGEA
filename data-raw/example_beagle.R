gl <- example_gl$likelihoods

fl <- strsplit(gl, ",")
fl1 <- unlist(lapply(fl, function(x) 10^((-1*as.numeric(x[1]))/10)))
fl2 <- unlist(lapply(fl, function(x) 10^((-1*as.numeric(x[2]))/10)))
fl3 <- unlist(lapply(fl, function(x) 10^((-1*as.numeric(x[3]))/10)))

fl <- cbind(fl1, fl2, fl3)
fl <- fl/rowSums(fl)



all(which(is.na(fl[,1])) == which(rowSums(is.na(fl)) == 3))
fl[is.na(fl)] <- 1/3


fl1 <- matrix(fl[,1], nrow = nrow(gl))
fl2 <- matrix(fl[,2], nrow = nrow(gl))
fl3 <- matrix(fl[,3], nrow = nrow(gl))


fl <- cbind(fl1, fl2, fl3)
k <- ncol(fl1)
fl <- fl[,order(rep(1:k, 3))]
fl <- as.data.table(fl)
cn <- rep(1:ncol(gl), each = 3)
cn <- paste0("Ind", formatC(cn, width = 3, format = "d", flag = "0"))
colnames(fl) <- cn
snpinfo <- data.table(marker = paste0("CHR1_", 1:nrow(fl)),
                      allele1 = 1,
                      allele2 = 2)

beagle <- cbind(snpinfo, fl)
data.table::fwrite(beagle, "inst/extdata/example.beagle.gz", sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE)
