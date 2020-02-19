#' @export
ploidy.inference <- function(x, chrom=NULL, start=NULL, end=NULL, penalty=25)
{
  # Make sure the penalties can be safely converted to a factor for splitting
  # purposes
  stopifnot(!any(duplicated(as.character(penalty))))
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(penalty))
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
  segments <-
    scquantum:::prof2invals(x, penalty, annotations, "chrom", "start", "end")
  if (is.null(chrom)) segments$chrom <- NULL
  filtered.segments <- segments[segments$length >= 20,]
  filtered.ratio.segments <- filtered.segments
  filtered.ratio.segments$mean <- filtered.ratio.segments$mean / mean(x)
  filtered.ratio.segments$se <- filtered.ratio.segments$se / mean(x)
  polar.quantogram <- Reduce(
    rbind,
    mapply(function(seg, l)
    {
      # Make sure that the order of the penalty isn't getting switched
      stopifnot(seg$penalty[1] == l)
      data.frame(
        penalty=l,
        s = seq(1, 8, length.out=100),
        polar_quantogram = scquantum:::weighted.ecf(
        seg$mean, seg$se,
        seq(1, 8, length.out=100))
      )
    },
    # Split while guaranteeing that the order stays the same
    split(filtered.ratio.segments,
          factor(as.character(filtered.ratio.segments$penalty),
                 levels=as.character(penalty))),
    penalty,
    SIMPLIFY=FALSE)
  )
  optimization.results <- Reduce(
    rbind,
    mapply(function(polar.quantogram, penalty)
    {
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
      return(data.frame(peak_location=peak.location,
                        peak_height=peak.height,
                        peak_phase=peak.phase))
    },
    split(polar.quantogram,
          factor(as.character(polar.quantogram$penalty), levels=as.character(penalty))),
    penalty,
    SIMPLIFY=FALSE)
  )
  # Construct the output list
  output <- with(optimization.results,
    list(penalty = penalty,
         multiply_ratios_by = peak_location,
         subtract_from_scaled_ratios = peak_phase / (2*pi),
         ploidy = peak_location - peak_phase / (2*pi),
         peak_height = peak_height,
         segmentation = segments,
         polar_quantogram = polar.quantogram,
         bincounts = bincounts)
  )
  class(output) <- c("scquantum_ploidy_inference", class(output))
  names(output$penalty) <- sprintf("penalty=%s", as.character(penalty))
  names(output$multiply_ratios_by) <- sprintf("penalty=%s", as.character(penalty))
  names(output$subtract_from_scaled_ratios) <- sprintf("penalty=%s", as.character(penalty))
  names(output$ploidy) <- sprintf("penalty=%s", as.character(penalty))
  names(output$peak_height) <- sprintf("penalty=%s", as.character(penalty))
  return(output)
}
