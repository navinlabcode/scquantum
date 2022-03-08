#' Ploidy inference with quantogram
#'
#' Infer ploidy of a cell, given a copy number profile.
#' Constructs a quantogram (either modular or cosine, depending on parameters).
#' The maximum of the quantogram is the estimated ploidy.
#' If unsegmented bincounts are given, segmentation will be performed using the
#' fused lasso.
#'
#' @param x Bincounts or segment means
#' @param chrom Optional chromosome numbers
#' @param start Optional bin/segment start positions
#' @param end Optional bin/segment end positions
#' @param penalty If segmenting, penalty parameter for the fused lasso (higher penalty, fewer segments)
#' @param do_segmentation Boolean, whether to do segmentation (set this to TRUE if giving unsegmented bincounts)
#' @param seg_length If giving already segmented data, length of each segment
#' @param iod If giving already segmented data, the index of dispersion of the bincount distribution (that is, within segments, not including between-segment variance)
#' @param mean_bincount If giving already segmented ratio values, the original mean bincount
#'
#' @return A ploidy inference object
#'
#' \describe{
#'   \item{penalty}{The segmentation penalty given as an argument, if any}
#'   \item{multiply_ratios_by}{To convert ratios to (unrounded) copy number estimates, multiply by this number}
#'   \item{subtract_from_scaled_ratios}{To convert ratios to (unrounded) copy numbers, after multiplying, subtract this number. Only required if the count data have some extra reads even at copy number 0, generally due to mapping problems}
#'   \item{ploidy}{The estimated ploidy}
#'   \item{peak_height}{The height of the quantogram peak at the estimated ploidy. Between 0 and 1. Higher values indicate a stronger signal}
#'   \item{segmentation}{The segmented values (either given as an argument, or produced interally by segmentation)}
#'   \item{polar_quantogram}{The complex-valued quantogram, whose absolute values measure consistency with each possible ploidy}
#'   \item{bincounts}{The raw bincounts given as an argument (if a segmentation was not given directly)}
#'   \item{theoretical_quantogram}{Based on the inferred copy numbers and index of dispersion, what the absolute value of the quantogram should look like. Deviation of this theoretical quantogram from the real one indicate that the ploidy estimate may be wrong}
#'   \item{theoretical_peak_height}{Height of the peak in the theoretical quantogram, measuring the expected strength of signal for the ploidy value}
#'   \item{confidence_ratio}{Ratio of actual to theoretical peak height. Values near (or above) 1 indicate the signal was as strong as would be expected gievn this data quality and ploidy; low values indicate that the ploidy inference may be wrong or that there are unexpected quality issues with the data}
#' }
#'
#' @examples
#' # Generating a random copy number profile
#' set.seed(705)
#' cns <- rpois(30, 3) + 1
#' x <- unlist(lapply(cns, function(cn) rpois(100, 25 * cn)))
#' annotations <- data.frame(chrom = 1, start = 1:length(x), end = 1:length(x))
#'
#' # Inferring ploidy
#' # Annotations and penalty are optional
#' estimate.from.bincounts <- ploidy.inference(x, annotations$chrom, annotations$start, penalty = 25)
#'
#' # Using scquantum internal functions to segment the data and estimate index
#' # of git@github.com:navinlabcode/scquantum.git
#' # dispersion
#' mu.est <- mean(x)
#' iod.est <- timeseries.iod(x)
#' seg <- prof2invals(x, 25, annotations, "chrom", "start", "end")
#' mean.est <- mean(x)
#' iod.est <- timeseries.iod(x)
#' estimate.from.segmentation <-
#'   ploidy.inference(
#'     seg$mean,
#'     seg$chrom,
#'     seg$start,
#'     seg$end,
#'     iod = iod.est,
#'     mean_bincount = mean.est,
#'     do_segmentation = FALSE
#'   )
#' @export
#'

ploidy.inference <- function(x, chrom = NULL, start = NULL, end = NULL, penalty = 25,
  do_segmentation = TRUE, seg_length = NULL, iod = NULL, mean_bincount = NULL)
{
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty) == 1)
  stopifnot(penalty > 0)
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
      prof2invals(x, penalty, annotations, "chrom", "start", "end")
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
    segments <- seg2invals(seg_mean, seg_length, iod, annotations)
  }
  if (is.null(chrom)) segments$chrom <- NULL
  filtered.segments <- segments[segments$length >= 20,]
  filtered.ratio.segments <- filtered.segments
  filtered.ratio.segments$mean <- filtered.ratio.segments$mean / mu.est
  filtered.ratio.segments$se <- filtered.ratio.segments$se / mu.est
  svals <- seq(1, 8, length.out=100)
  polar.quantogram <- data.frame(
        s = svals,
        polar_quantogram = weighted.ecf(
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
    theoretical_quantogram = expected.peak.heights(
      cn.est,
      filtered.ratio.segments$se,
      optimization.results$peak_location,
      svals))
  theoretical.peak.height <- expected.peak.heights(
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
