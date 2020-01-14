### Functions for processing raw profiles

# A function for estimating the index of dispersion, which is used when
# estimating standard errors for each segment mean. This is an old version,
# which I'm keeping just to compare to the newer version, "timeseries.iod",
# below.

#' @export
diffsum.tm.2 <- function(v)
{
  y <- v[-1]
  x <- v[-length(v)]
  vals <- (y-x)^2/(y+x+1)
  q <- quantile(vals, probs=.95)
  tm <- mean(vals[vals <= q])
  mean.est <- 1.05*(.95*tm + .05 * q)
  return(mean.est)
}

# Estimate the standard deviation of a normal distribution with L2E, assuming
# the normal distribution has mean 0.

#' @export
l2e.normal.sd <- function(xs)
{
  optim.result <- optimize(
    # L2E loss function
    f=function(sd)
    # "Data part", the sample average of the likelihood
    -2 * mean(dnorm(xs, sd=sd)) +
    # "Theta part", the integral of the squared density
      1/(2*sqrt(pi)*sd),
    # Parameter: standard deviation of the normal distribution fit
    interval = c(0, diff(range(xs))))
  return(optim.result$minimum)
}

# A function for estimating the index of dispersion, which is used when
# estimating standard errors for each segment mean

#' @export
timeseries.iod <- function(v)
{
  # Differences between pairs of values
  y <- v[-1]
  x <- v[-length(v)]
  # Normalize the differences using the sum. The result should be around zero,
  # plus or minus square root of the index of dispersion
  vals.unfiltered <- (y-x)/sqrt(y+x)
  vals <- vals.unfiltered[y + x  >= 1]
  # Assuming most of the normalized differences follow a normal distribution,
  # estimate the standard deviation
  val.sd <- l2e.normal.sd(vals)
  # Square this standard deviation to obtain an estimate of the index of
  # dispersion
  return(val.sd^2)
}

# A variance-stabilizing transform for data where the mean varies but the index
# of dispersion stays the same. "gat" stands for "generalized Anscombe transform"

#' @export
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
tf.dp <- Vectorize(function(y, lam)
{
  n <- length(y)
  stopifnot(length(lam)==1)
  .C("tf_dp_wrapper", n=as.integer(n), y=as.numeric(y), lam=as.double(lam)[1], beta=numeric(n))$beta
}, "lam")

# Segment a profile, and summarize the segments. Output has columns lambda,
# mean, se, length, start_index, and end_index. Segmentation is performed for
# each input penalty value ("lambda"), and each row int he output represents a
# segment, giving the estimated mean of the segment, the standard error of the
# mean, and the number of elements which were in the segment.
segment.summarize <- function(inprof, lambda, trans, seg, loc, se)
{
  # If the input penalties don't have names, name them. This is intended to make
  # grouping by penalty safe without worrying about flaoting point problems.
  # Transform to try and get data which are normally distributed around a
  # segment-specific mean, with the same variance for each segment
  transformed.profile <- trans(inprof)

  # For each penalty value, get a segmented profile in long format--that is, of
  # the same length of the original profile, but piecewise constant. Each
  # segmented profile, corresponding to a penalty value, will be a column of a
  # matrix.
  nlambda <- length(lambda)
  if (nlambda > 1)
  {
    segmented.profiles <- seg(transformed.profile, lambda)
  } else if (nlambda == 1)
  # If lambda is of length 1, by default we will get an array, not a matrix, and
  # it needs to be converted to a matrix for downstream stuff to work
  {
    segmented.profiles <- as.matrix(seg(transformed.profile, lambda), ncol=1)
  } else
  # What's the remaining case? Length is 0. Conceptually I could return a table
  # with no rows in that case, but the function is not really supposed to do
  # nothing and not segment anything, so I'll make that an error case
  {
    stop("No penalty values given")
  }

  # For each penalty value, number elements of the profile according to what
  # segment they're in
  segnums <- apply(segmented.profiles, 2, function(segmented.profile)
    cumsum(c(TRUE, abs(diff(segmented.profile)) > 0.1))
  )

  # Summarize the segments, recording three numbers: estimate the mean of each
  # segment, the standard error of the mean, and record the length of the
  # segment.
  means <- lapply(1:nlambda, function(i)
    tapply(inprof, segnums[,i], loc)
  )
  standard.errors <- lapply(1:nlambda, function(i)
    tapply(inprof, segnums[,i], se)
  )
  lengths <- lapply(1:nlambda, function(i)
    tapply(inprof, segnums[,i], length)
  )

  # From the lengths, get the start and end indices
  end.indices <- lapply(lengths, cumsum)
  start.indices <- lapply(end.indices, function(x) c(0, x[-length(x)]) + 1)

  Reduce(rbind, lapply(1:length(lambda), function(i)
    data.frame(penalty = lambda[i], mean=means[[i]], se=standard.errors[[i]],
               length=lengths[[i]], start_index = start.indices[[i]], end_index = end.indices[[i]])
  ))
}

#' @export
prof2invals <- function(
  # The input required for segmentation: the profile to be segmented, and the
  # penalty values for the segmentation
  inprof, lambda,
  # The input required for annotation: the data frame containing the
  # annotations, and the names of the columns which are going to be used. This
  # also affects the segmentation, since chromosome boundaries are always
  # segment boundaries.
  annotations, chrom.colname, bin.start.colname, bin.end.colname)
{
  # Estimate the index of dispersion, which will be used for estimating standard
  # errors of segment means, and for transformation
  iod.est <- scquantum:::timeseries.iod(inprof)

  # Split the annotations into two parts, one which will be accessed using the
  # segment start indices, and the other which will be accessed using the
  # segment end indices. The chromosome column can go in either one of these,
  # but of course the bin start positions need to go in the first, and the bin
  # end positions need to go in the second.
  left.annotations <- annotations[,c(chrom.colname, bin.start.colname),drop=FALSE]
  right.annotations <- annotations[,bin.end.colname,drop=FALSE]

  # Segment the profile and annotate the results
  annotated.segmented.counts <- Reduce(rbind, mapply(function(segments, annotations1, annotations2)
    cbind(annotations1[segments$start_index,,drop=FALSE], annotations2[segments$end_index,,drop=FALSE], segments),
    # Split the profile by chromosome and segment each chromosome separately,
    # guaranteeing that chromosome boundaries are segment boundaries
    tapply(inprof, annotations[[chrom.colname]], scquantum:::segment.summarize, lambda=lambda,
      # The functions to be used to transform the data, to segment it, to estimate
      # the segment means, and to estimate the standard error of the means
      trans=function(x) log2(x / mean(x) + 0.15),
      seg = scquantum:::tf.dp, loc = median,
      se = function(x) sqrt(pi/2) * sqrt(iod.est * median(x) / length(x))
    ),
    # Chromosome annotations and bin start annotations, split by chromosome
    split(left.annotations, annotations[[chrom.colname]]),
    # Bin end annotations, split by chromosome
    split(right.annotations, annotations[[chrom.colname]]),
    # Can be read either as a statement that the output is to be returned as a
    # list, or a sarcastic comment about how the code was written
    SIMPLIFY=FALSE))

  return(annotated.segmented.counts[,c("penalty", chrom.colname, bin.start.colname, bin.end.colname, "mean", "se", "length")])
}

### Empirical characteristic functions and maxima

# Evaluate the weighted empirical characteristic function
#' @export
weighted.ecf <- Vectorize(function(y, sds, s)
{
  variances <- 1 - exp(-4 * pi^2 * sds^2 * s^2)
  weights <- (1/variances) / sum(1/variances)
  sum(weights * exp(1i * (2*pi) * s * y))
}, 's')

# Global maximum of the weighted empirical characteristic function
#' @export
ecf.global.max <- function(y, sds, smin=1, smax = 8)
{
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

# Local maximum of the weighted empirical characteristic function
#dyn.load("ecf_local_max.dylib")
ecf.local.max <- function(s, y, sds)
{
  stopifnot(length(y) == length(sds))
  if (any(is.na(c(s,y,sds)) | is.nan(c(s,y,sds)) | is.infinite(c(s, y, sds))))
  {
    return(c(peak_location=NA, peak_height=NA, peak_phase=NA))
  }
  call.results <-
     .C("ecf_local_max", s=as.double(s), nys=length(y), ys=as.double(y),
     nsigmas=length(sds), sigmas=as.double(sds),
     peak_location=double(1), peak_height=double(1), peak_phase=double(1))
  return(unlist(call.results[c("peak_location", "peak_height", "peak_phase")]))
}
