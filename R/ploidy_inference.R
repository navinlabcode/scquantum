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