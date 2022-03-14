#' Functions for processing raw profiles




# Estimate the standard deviation of a normal distribution with L2E, assuming
# the normal distribution has mean 0.


#' @export
#' @keywords internal
#' @name l2e.normal.sd
#' @docType methods
#' @rdname internals
l2e.normal.sd <- function(xs)
{
  # Need at least two values to get a standard deviation
  stopifnot(length(xs) >= 2)
  optim.result <- stats::optimize(
    # L2E loss function
    f=function(sd)
    # "Data part", the sample average of the likelihood
    -2 * mean(stats::dnorm(xs, sd=sd)) +
    # "Theta part", the integral of the squared density
      1/(2*sqrt(pi)*sd),
    # Parameter: standard deviation of the normal distribution fit
    interval = c(0, diff(range(xs))))
  return(optim.result$minimum)
}


# A variance-stabilizing transform for data where the mean varies but the index
# of dispersion stays the same. "gat" stands for "generalized Anscombe transform"

#' @export
#' @keywords internal
#' @name gat
#' @docType methods
#' @rdname internals
gat <- function(x, iod)
{
  # When index of dispersion is above 4, transformation will be imaginary for
  # low inputs. So, you need to check if it's real-valued. If it's imaginary,
  # return NA.
  ifelse(x + (4-iod)/8 < 0,
    NA,
    # Due to vectorization, this will get executed even when the result is
    # imaginary, resulting in a warning. Therefore, suppress warnings.
    suppressWarnings((2 / sqrt(iod)) * sqrt(x + (4-iod)/8))
  )
}

#' @useDynLib scquantum tf_dp_wrapper
#' @export
#' @keywords internal
#' @name tf.dp
#' @docType methods
#' @rdname internals
tf.dp <- Vectorize(function(y, lam)
{
  n <- length(y)
  stopifnot(length(lam)==1)
  .C("tf_dp_wrapper", n=as.integer(n), y=as.numeric(y),
     lam=as.double(lam)[1], beta=numeric(n))$beta
}, "lam")

# Segment a profile, and summarize the segments. Output has columns penalty,
# mean, se, length, start_index, and end_index. Each row in the output
# represents a segment, giving the estimated mean of the segment, the standard
# error of the mean, and the number of elements which were in the segment.
#' @export
#' @keywords internal
#' @name segment.summarize
#' @docType methods
#' @rdname internals
segment.summarize <- function(inprof, penalty, trans, seg, loc, se)
{
  # If the input penalties don't have names, name them. This is intended to make
  # grouping by penalty safe without worrying about flaoting point problems.
  # Transform to try and get data which are normally distributed around a
  # segment-specific mean, with the same variance for each segment
  transformed.profile <- trans(inprof)
  stopifnot(all(!is.nan(transformed.profile)) & all(!is.na(transformed.profile)))

  segmented.profile <- seg(transformed.profile, penalty)

  # For each penalty value, number elements of the profile according to what
  # segment they're in
  segnums <- cumsum(c(TRUE, abs(diff(segmented.profile)) > 0.1))

  # Summarize the segments, recording three numbers: estimate the mean of each
  # segment, the standard error of the mean, and record the length of the
  # segment. Values are not transformed, since they are obtained from the
  # original profile, using only the segment numbering derived from the
  # segmented transformed profile.
  means <- tapply(inprof, segnums, loc)
  standard.errors <- tapply(inprof, segnums, se)
  lengths <- tapply(inprof, segnums, length)

  # From the lengths, get the start and end indices
  end.indices <- cumsum(lengths)
  start.indices <- c(0, end.indices[-length(end.indices)]) + 1

  data.frame(mean=means, se=standard.errors,
             length=lengths, start_index = start.indices,
             end_index = end.indices)
}

#' @export
#' @keywords internal
#' @name prof2invals
#' @docType methods
#' @rdname internals
prof2invals <- function(
  # The input required for segmentation: the profile to be segmented, and the
  # penalty value for the segmentation
  inprof, penalty,
  # The input required for annotation: the data frame containing the
  # annotations, and the names of the columns which are going to be used. This
  # also affects the segmentation, since chromosome boundaries are always
  # segment boundaries.
  annotations, chrom.colname, bin.start.colname, bin.end.colname)
{

  stopifnot(typeof(chrom.colname) == "character")
  stopifnot(typeof(bin.start.colname) == "character")
  stopifnot(typeof(bin.end.colname) == "character")
  stopifnot("data.frame" %in% class(annotations))

  # Estimate the index of dispersion, which will be used for estimating standard
  # errors of segment means, and for transformation
  iod.est <- timeseries.iod(inprof)

  # Split the annotations into two parts, one which will be accessed using the
  # segment start indices, and the other which will be accessed using the
  # segment end indices. The chromosome column can go in either one of these,
  # but of course the bin start positions need to go in the first, and the bin
  # end positions need to go in the second.
  left.annotations <- annotations[,c(chrom.colname, bin.start.colname),drop=FALSE]
  # If you don't want to annotate with bin ends, you can put in NULL, and the
  # right annotations will be a data frame with no columns
  right.annotations <- if (is.null(bin.end.colname))
  {
    data.frame()[1:nrow(annotations),]
  } else
  {
      out <- annotations[,bin.end.colname,drop=FALSE]
      # When bin.start.colname overlaps with bin.end.colname, add a suffix to
      # the bin end column names to avoid duplicate column names
      while (any(colnames(out) %in% colnames(left.annotations)))
      {
        colnames(out)[colnames(out) %in% colnames(left.annotations)] <-
          sapply(colnames(out)[colnames(out) %in% colnames(left.annotations)],
                 paste, "end", sep=".")
      }
      out
  }
  # Segment the profile and annotate the results
  annotated.segmented.counts <-
    Reduce(rbind, mapply(function(segments, annotations1, annotations2)
    cbind(annotations1[segments$start_index,,drop=FALSE],
          annotations2[segments$end_index,,drop=FALSE],
          segments),
    # Split the profile by chromosome and segment each chromosome separately,
    # guaranteeing that chromosome boundaries are segment boundaries
    tapply(inprof, annotations[[chrom.colname]], segment.summarize,
           penalty=penalty,
      # The functions to be used to transform the data, to segment it, to estimate
      # the segment means, and to estimate the standard error of the means
      trans = function(x)
      {
        gat.result <- gat(x, iod=iod.est);
        ifelse(is.na(gat.result), 0, gat.result)
      },
      seg = tf.dp, loc = stats::median,
      se = function(x) sqrt(pi/2) * sqrt(iod.est * stats::median(x) / length(x))
    ),
    # Chromosome annotations and bin start annotations, split by chromosome
    split(left.annotations, annotations[[chrom.colname]]),
    # Bin end annotations, split by chromosome
    split(right.annotations, annotations[[chrom.colname]]),
    # Can be read either as a statement that the output is to be returned as a
    # list, or a sarcastic comment about how the code was written
    SIMPLIFY=FALSE))

  return(annotated.segmented.counts[,c(
    colnames(left.annotations), colnames(right.annotations),
    "mean", "se", "length"
  )])
}

#' @export
#' @keywords internal
#' @name seg2invals
#' @docType methods
#' @rdname internals
seg2invals <- function(seg_mean, seg_length, iod, annotations)
{
  se <- sqrt(iod * seg_mean / seg_length)
  return(cbind(annotations, data.frame(mean = seg_mean, se = se, length = seg_length)))
}


### Empirical characteristic functions and maxima

# Evaluate the weighted empirical characteristic function
#' @export
#' @keywords internal
#' @name weighted.ecf
#' @docType methods
#' @rdname internals
weighted.ecf <- Vectorize(function(y, sds, s)
{
  stopifnot(length(y) == length(sds))
  y <- y[sds > 0.00001]
  sds <- sds[sds > 0.00001]
  means <- exp(-2 * pi^2 * sds^2 * s^2)
  variances <- 1 - exp(-4 * pi^2 * sds^2 * s^2)
  unnormalized.weights <- means / variances
  weights <- unnormalized.weights / sum(unnormalized.weights)
  sum(weights * exp(1i * (2*pi) * s * y))
}, 's')

# Global maximum of the weighted empirical characteristic function
#' @export
#' @keywords internal
#' @name ecf.global.max
#' @docType methods
#' @rdname internals
ecf.global.max <- function(y, sds, smin=1, smax = 8)
{
  stopifnot(length(y) == length(sds))
  svals <- seq(smin, smax, by=0.01)
  polar.quantogram <- weighted.ecf(y, sds, svals)
  max.index <- which.max(Mod(polar.quantogram))
  if (length(max.index) == 1)
  {
    return(c(
      peak_location = svals[max.index],
      peak_height = Mod(polar.quantogram)[max.index],
      peak_phase = Arg(polar.quantogram)[max.index]
    ))
  } else
  {
    return(c(
      peak_location = NA,
      peak_height = NA,
      peak_phase = NA
    ))
  }
}

# Takes in copy numbers, ratio standard errors, the scaling factor, and input
# values to the characteristic function. Outputs absolute value of the
# characteristic function
#' @export
#' @keywords internal
#' @name expected.peak.heights
#' @docType methods
#' @rdname internals
expected.peak.heights <- function(cn, ratio.se, scale, svals)
{
  nonzero.cn <- cn[cn > 0]
  nonzero.se <- (ratio.se * scale)[cn > 0]

  # Produce matrices, where each row is an s value, and each column is a segment
  # Characteristic function of the noise
  noise.cf <- t(sapply(svals, function(s)
    exp(-(1/2) * nonzero.se^2 * (2 * pi * s)^2)
  ))
  # Characteristic function of the copy number. I'm considering the copy number
  # to be a constant, but it still has a characteristic function, although a
  # trivial one
  # I'm leaving out the offset for now because I'm just calculating the modulus,
  # which is not affected by it; but I should put it in just so the definitions
  # don't look wrong
  cn.cf <- t(sapply(svals, function(s)
    exp(1i * (nonzero.cn / scale) * (2 * pi * s))
  ))
  # E[X] = Psi(2 pi s)
  expected.values <- cn.cf * noise.cf
  mean.norms <- noise.cf
  # Var(X) = 1 - | Psi(2 pi s) |^2
  variances <- 1 - mean.norms^2
  unnormalized.weights <- mean.norms / variances
  partition.function <- apply(unnormalized.weights, 1, sum)
  weights <- unnormalized.weights / partition.function
  weights[apply(weights, 1, function(x) any(is.nan(x)))] <- NA
  this.theoretical.quantogram <-
    apply(weights * expected.values, 1, sum)
  background <- apply(weights^2 * variances, 1, sum)
  expected.peak.height.function <- sqrt(Mod(this.theoretical.quantogram)^2 + background)
  return(expected.peak.height.function)
}

#' @export
#' @keywords internal
#' @name quantum.sd
#' @docType methods
#' @rdname internals
quantum.sd <- function(x, mu)
{
  y <- x / mu
  svals <- seq(1, 8, length.out=700)
  polar.quantogram <- sapply(svals, function(s) mean(exp(1i * (2*pi) * s * y)))
  peak.height <- Mod(polar.quantogram)[which.max(Mod(polar.quantogram))[1]]
  peak.location <- svals[which.max(Mod(polar.quantogram))[1]]
  ratio.variance <- (log(peak.height) / (-2 * pi^2)) * peak.location^2
  ratio.variance * mu
}

#' @export
#' @keywords internal
#' @name irises.pluspurple
#' @docType methods
#' @rdname internals
irises.pluspurple <-
  c(`0`="#25292E", `1`="#9F3D0C", `2`="#CF601F", `3`="#E1DE96",
    `4`="#89AD71", `5`="#89C9E0", `6`="#556bb7", `>=7`="#8b44bb")

# A function for estimating the index of dispersion, which is used when
# estimating standard errors for each segment mean
#' @export
#' @keywords internal
#' @name timeseries.iod
#' @docType methods
#' @rdname internals
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
