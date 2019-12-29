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
    lower=mad(xs)/100, upper=mad(xs)*100)
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

# Version of prof2invals which records more information
# Output table has chrom, seg_start_relpos, seg_end_relpos, length, location,
# ratio, and ratio_sd
prof2invals <- function(inprof, lambda)
{  
  # Transform to try and get data which are normally distributed around a
  # segment-specific mean, with the same variance for each segment
  transformed.profile <- log2(inprof/mean(inprof) + 0.15)
  # Split up on chromosome boundaries (it should really be arm boundaries, but
  # that can wait)
  chrom.profiles <- split(transformed.profile, chrom.column)
  # For each penalty value, get a segmented profile in long format
  segmented.profiles <-
    Reduce(cbind, lapply(chrom.profiles, function(x) 
                         flsaGetSolution(flsa(x), lambda2=lambda)))
  # For each penalty value, number bins according to what segment they're in
  is.newseg <- t(apply(segmented.profiles, 1, function(segmented.profile)
                       c(TRUE, (abs(diff(segmented.profile)) > 0.1) | (diff(chrom.column) != 0))))
  segnum <- t(apply(is.newseg, 1, cumsum))
  # Get a location estimate for each segment, and the length of each segment
  rownames(segnum) <- sprintf("lambda%.5f", lambda)
  loc.length.table <- as_tibble(t(segnum)) %>%
    mutate(bincounts=inprof) %>%
    bind_cols(bin.start.end) %>%
    gather(key="penalty", value="segnum", -bincounts,
           -chrom, -bin.start, -bin.end, -bin.number,
           -abs_start, -abs_end) %>%
    group_by(penalty, segnum) %>%
    summarize(length=length(bincounts), location=median(bincounts),
              start_chrom = min(chrom), end_chrom = max(chrom),
              start_pos = min(bin.start), end_pos = max(bin.end),
              abs_start = min(abs_start), abs_end = max(abs_end)) %>%
    ungroup()
  # Calculate index of dispersion of the profile
  iod <- diffsum.tm.2(inprof)
  # Calculate mean of the part of the profile that's included
  region.total <- with(loc.length.table,
                       sum(length * location) / sum(length))
  # Ratios and standard deviation of the ratios
  loc.length.weight.ratio.table <-
    dplyr::mutate(loc.length.table,
      ratio=location / region.total,
      ratio_sd = sqrt(pi/2) * sqrt(iod * location / length) / region.total)
  return(loc.length.weight.ratio.table)
}

### Empirical characteristic functions and maxima

# Evaluate the weighted empirical characteristic function
weighted.ecf <- Vectorize(function(y, sds, s)
{
  variances <- 1 - exp(-4 * pi^2 * sds^2 * s^2)
  weights <- (1/variances) / sum(1/variances)
  sum(weights * exp(1i * (2*pi) * s * y))
}, 's')

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
