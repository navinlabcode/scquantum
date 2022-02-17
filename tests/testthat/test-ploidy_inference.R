test_that("Runs without error", {
  set.seed(705)
  cns <- rpois(30, 3) + 1
  x <- unlist(lapply(cns, function(cn) rpois(100, 25 * cn)))
  # Inferring ploidy
  estimate.from.bincounts <- ploidy.inference(x)
  expect_length(estimate.from.bincounts$peak_height, 1)
  annotations <- data.frame(chrom = 1, start = 1:length(x), end = 1:length(x))
  # Using scquantum internal functions to segment the data and estimate index of
  # dispersion
  mu.est <- mean(x)
  iod.est <- scquantum:::timeseries.iod(x)
  seg <- scquantum:::prof2invals(x, 25, annotations, "chrom", "start", "end")
  mean.est <- mean(x)
  iod.est <- scquantum:::timeseries.iod(x)
  estimate.from.segmentation <- ploidy.inference(seg$mean, seg$chrom, seg$start, seg$end, iod = iod.est, mean_bincount = mean.est, do_segmentation = FALSE)
  expect_length(estimate.from.segmentation$peak_height, 1)
})
