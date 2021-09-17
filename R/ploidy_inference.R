#' @export
ploidy.inference <- function(x, chrom = NULL, start = NULL, end = NULL, penalty = 25,
  do_segmentation = TRUE, seg_length = NULL, iod = NULL)
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
    segments <-
      scquantum:::prof2invals(x, penalty, annotations, "chrom", "start", "end")
  } else
  {
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
  filtered.ratio.segments$mean <- filtered.ratio.segments$mean / mean(x)
  filtered.ratio.segments$se <- filtered.ratio.segments$se / mean(x)
  polar.quantogram <- data.frame(
        penalty=penalty,
        s = seq(1, 8, length.out=100),
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
  return(output)
}

#' @export
plot.scquantum_ploidy_inference <- function(ploidy.inference)
{
  # Segmented profile
  mean.bincount <- mean(ploidy.inference$bincounts$bincount)
  bincount2cn <- function(x)
  {
    with(ploidy.inference,
         (x / mean.bincount) * multiply_ratios_by - subtract_from_scaled_ratios)
  }
  profile.plotelems <- if (length(ploidy.inference$penalty) == 1)
  {
    cn.string <- with(ploidy.inference, ifelse(
      round(bincount2cn(bincounts$bincount)) >= 7,
      ">=7",
      as.character(round(bincount2cn(bincounts$bincount)))
    ))
    cn.factor <- factor(cn.string, levels=c(as.character(0:6), ">=7"))
    bincounts.with.cn.estimates <- ploidy.inference$bincounts
    bincounts.with.cn.estimates$cn <- cn.factor
    segmentation.for.plot <- ploidy.inference$segmentation
    segmentation.for.plot <- segmentation.for.plot[segmentation.for.plot$length >= 20,]
    list(
      geom_point(aes(x=pos, y=bincount, colour=cn),
                 data=bincounts.with.cn.estimates,
                 alpha=0.5, size=0),
      geom_segment(aes(x=start, xend=end, y=mean, yend=mean),
                   data=segmentation.for.plot),
    scale_colour_manual(
      values = scquantum:::irises.pluspurple,
      guide=guide_legend(title="copy\nnumber", override.aes=list(alpha=1, size=2))),
    scale_y_continuous(sec.axis = sec_axis(bincount2cn, name="copy number", breaks=0:10),
                       limits=c(0, max(ploidy.inference$segmentation$mean[ploidy.inference$segmentation$length >= 20]) * 1.1))
    )
  } else
  {
    list(
      geom_point(aes(x=pos, y=bincount),
                 data=ploidy.inference$bincounts, alpha=0.1, size=0),
      geom_segment(aes(x=start, xend=end, y=mean, yend=mean, colour=as.factor(penalty)),
                   data=ploidy.inference$segmentation),
      scale_colour_discrete(guide=guide_legend(title="penalty")),
      ylim(0, max(ploidy.inference$segmentation$mean[ploidy.inference$segmentation$length >= 20]) * 1.1)
    )
  }
  if (!is.null(ploidy.inference$segmentation$chrom))
  {
    profile.plotelems <- c(profile.plotelems,
      facet_grid(cols=vars(factor(chrom, levels=stringr::str_sort(unique(chrom), numeric=TRUE))),
                 space="free_x", scales="free_x"))
  }
  segmented.profile <- Reduce(`+`, c(list(ggplot(ploidy.inference$bincounts)), profile.plotelems)) +
    cowplot::theme_cowplot() +
    theme(panel.spacing.x=unit(0, "in"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text=element_text(angle=90, size=6)) +
    xlab("position")

  # Histogram or frequency polygon
  binwidth.rule <- function(x)
  {
    default.sd <- mean.bincount / 100
    est.sd <- scquantum:::quantum.sd(x, mean.bincount)
    # Not the right length
    return(3.5 * max(default.sd, est.sd) / (length(x)/5)^(1/3))
  }
  distribution.geom <-
    if (length(ploidy.inference$penalty)==1)
    {
      geom_histogram(
        aes(mean, stat(count)),
        binwidth=3.5 * median(ploidy.inference$segmentation$se[ploidy.inference$segmentation$length >= 20]) /
          (sum(ploidy.inference$segmentation$length >= 20)/ploidy.inference$ploidy)^(1/3)
      )
    } else
    {
      geom_freqpoly(
        aes(mean, stat(count), colour=as.factor(penalty)),
        binwidth=binwidth.rule
      )
    }
  distribution <-
    ggplot(ploidy.inference$segmentation[ploidy.inference$segmentation$length >= 20,]) +
    cowplot::theme_cowplot() +
    distribution.geom +
    theme(legend.position='none') +
    xlab("segment mean") + ylab("number of segments")

  # Modular quantogram
  quantogram.geom <- if (length(ploidy.inference$penalty) == 1)
  {
    geom_line(aes(x=1 / s * mean.bincount,
                  y=Mod(polar_quantogram)))
  } else
  {
    geom_line(aes(x=1 / s * mean.bincount,
                  y=Mod(polar_quantogram), colour=as.factor(penalty)))
  }
  quantogram <- ggplot(ploidy.inference$polar_quantogram) +
    cowplot::theme_cowplot() +
    quantogram.geom +
    theme(legend.position='none') +
    xlab("reads per copy") +
    scale_y_continuous(name="score", limits=c(0,1))

  return(segmented.profile / (distribution | quantogram))
}

# A function for estimating the index of dispersion, which is used when
# estimating standard errors for each segment mean

#' @export
timeseries.iod <- function(v)
{
  # 3 elements, 2 differences, can find a standard deviation
  stopifnot(length(v) >= 3)
  # Differences between pairs of values
  y <- v[-1]
  x <- v[-length(v)]
  # Normalize the differences using the sum. The result should be around zero,
  # plus or minus square root of the index of dispersion
  vals.unfiltered <- (y-x)/sqrt(y+x)
  # Remove divide by zero cases, and--considering this is supposed to be count
  # data--divide by almost-zero cases
  vals <- vals.unfiltered[y + x  >= 1]
  # Check that there's anything left
  stopifnot(length(vals) >= 2)
  # Assuming most of the normalized differences follow a normal distribution,
  # estimate the standard deviation
  val.sd <- l2e.normal.sd(vals)
  # Square this standard deviation to obtain an estimate of the index of
  # dispersion
  return(val.sd^2)
}
