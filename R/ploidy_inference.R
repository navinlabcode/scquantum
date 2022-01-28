#' @export
ploidy.inference <- function(x, chrom = NULL, start = NULL, end = NULL, penalty = 25,
  do_segmentation = TRUE, seg_length = NULL, iod = NULL, mean_bincount = NULL)
{
  # Make sure the penalties can be safely converted to a factor for splitting
  # purposes
  stopifnot(!any(duplicated(as.character(penalty))))
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty) == 1)
  if (!is.null(chrom))
  {
    annotations <- data.frame(chrom=chrom)
  } else
  {
    annotations <- data.frame(chrom=rep.int("dummy_chrom", length(x)))
  }
  if (!is.null(start))
  {
    annotations$start <- start
    if (!is.null(end))
    {
      annotations$end <- end
    } else
    {
      # Use bin start values as bin end values if no end values are given
      annotations$end <- start
    }
  } else
  {
    annotations$start <- 1:length(x)
    annotations$end <- 1:length(x)
  }
  bincounts <- data.frame(bincount=x)
  if (!is.null(chrom)) bincounts$chrom <- chrom
  if (!is.null(start))
  {
    bincounts$pos <- start
  } else
  {
    bincounts$pos <- 1:length(x)
  }
  if (do_segmentation)
  {
    mu.est <- mean(x)
    segments <-
      scquantum:::prof2invals(x, penalty, annotations, "chrom", "start", "end")
  } else
  {
    mu.est <- mean_bincount
    seg_mean <- x
    stopifnot(!is.null(iod))
    stopifnot(!is.null(seg_length) | (!is.null(start) & !is.null(end)))
    if (is.null(seg_length))
    {
      seg_length <- end - start + 1
    }
    segments <- scquantum:::seg2invals(seg_mean, seg_length, iod, annotations)
  }
  if (is.null(chrom)) segments$chrom <- NULL
  filtered.segments <- segments[segments$length >= 20,]
  filtered.ratio.segments <- filtered.segments
  filtered.ratio.segments$mean <- filtered.ratio.segments$mean / mu.est
  filtered.ratio.segments$se <- filtered.ratio.segments$se / mu.est
  svals <- seq(1, 8, length.out=100)
  polar.quantogram <- data.frame(
        s = svals,
        polar_quantogram = scquantum:::weighted.ecf(
        filtered.ratio.segments$mean, filtered.ratio.segments$se,
        seq(1, 8, length.out=100))
      )
  optimization.results <- {
      max.index <- which.max(Mod(polar.quantogram$polar_quantogram))
      if (length(max.index) == 1) {
          peak.location <- polar.quantogram$s[max.index]
          peak.height <- Mod(polar.quantogram$polar_quantogram)[max.index]
          peak.phase <- Arg(polar.quantogram$polar_quantogram)[max.index]
      } else
      {
          peak.location <- NA
          peak.height <- NA
          peak.phase <- NA
      }
      data.frame(peak_location=peak.location,
                 peak_height=peak.height,
                 peak_phase=peak.phase)
    }
  cn.est <- round(filtered.ratio.segments$mean * optimization.results$peak_location - optimization.results$peak_phase / (2*pi))
  theoretical.quantogram <- data.frame(s = svals,
    theoretical_quantogram = scquantum:::expected.peak.heights(
      cn.est,
      filtered.ratio.segments$se,
      optimization.results$peak_location,
      svals))
  theoretical.peak.height <- scquantum:::expected.peak.heights(
      cn.est,
      filtered.ratio.segments$se,
      optimization.results$peak_location,
      optimization.results$peak_location)
  # Construct the output list
  output <- with(optimization.results,
    list(penalty = penalty,
         multiply_ratios_by = peak_location,
         subtract_from_scaled_ratios = peak_phase / (2*pi),
         ploidy = peak_location - peak_phase / (2*pi),
         peak_height = peak_height,
         segmentation = segments,
         polar_quantogram = polar.quantogram,
         bincounts = bincounts,
         theoretical_quantogram = theoretical.quantogram,
         theoretical_peak_height = theoretical.peak.height,
         confidence_ratio = peak_height / theoretical.peak.height)
  )
  class(output) <- c("scquantum_ploidy_inference", class(output))
  return(output)
}
