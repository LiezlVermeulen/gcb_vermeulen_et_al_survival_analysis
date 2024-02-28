## ----------------------------------------------------
## Script for detecting and characterising breakpoints
## ----------------------------------------------------
# Load required packages
library(rgdal) 
library(RColorBrewer)
library(raster)
library(bfast)
library(sf)
library(signal)
library(zoo)
library(ggplot2)

wd <- 'C:/Users/u0142455/Documents/PhD/Processing/ch1/SA'
setwd(wd)

# Load datasets
aoi <- shapefile("C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/aoi/savanna_studyarea.shp") # study area
rue<-brick("data/rue/RUE_NVDI3gLIN_CHIRPSoptisum.grd") # RUE time series
rue_sc <- calc(rue,fun = function(x) {
  x/10000 # rescale according to scale factor
} )
plot(rue_sc)

### Functions
## Applying bfast01classify to a whole raster
bf01c_raster<-function(raster,start,end,frequency,formula=response~trend+harmon,thrs=nlayers(raster)/2,h=0.15){
  r<-raster
  m<-as.matrix(r)
  ft=c()
  fst=c()
  fs=c()
  for (x in 1:ncell(r)){
    if(sum(is.na(m[x,]))<thrs){  
      pix=c()
      for (i in 1:nlayers(r)){
        pix[i]<-as.vector(m[x,i])
      }
      ts<-ts(pix,start=start,end=end,frequency=frequency)
      b1<-bfast01(ts, formula = NULL, test = c("OLS-MOSUM", "BIC"), level = 0.05, aggregate = any,
                  trim = 0.2, bandwidth = h, functional = "max",
                  order = 0, lag = NULL, slag = NULL, na.action = na.omit, stl = "none")
      b1c<-bfast01classify(b1, typology="drylands")
      ft[x]<-b1c$flag_type
      fst[x]<-b1c$flag_subtype
      fs[x]<-b1c$flag_significance
    }else{
      ft[x]<-NA
      fst[x]<-NA
      fs[x]<-NA
    }
  }
  results<-list(ft,fst,fs)
  names(results)<-c("flag_type","flag_subtype","flag_sign")
  
  type_raster<-raster(extent(r),nrow=nrow(r),ncol=ncol(r),crs=proj4string(r))
  values(type_raster)<-results[[1]]
  
  subtype_raster<-raster(extent(r),nrow=nrow(r),ncol=ncol(r),crs=proj4string(r))
  values(subtype_raster)<-results[[2]]
  
  sign_raster<-raster(extent(r),nrow=nrow(r),ncol=ncol(r),crs=proj4string(r))
  values(sign_raster)<-results[[3]]
  
  output<-list(type_raster,subtype_raster,sign_raster)
  return(output)
}

## Running bfast01 with different bandwidths
for(i in 2:33) {
  bf<-bf01c_raster(rue_crop,start=c(1982),end=c(2015),frequency=1,formula=response~trend,thrs=12,h=i/34)
  #out<-stack(bf[[1]],bf[[2]],bf[[3]])
  type <- bf[[1]]
  subtype <- bf[[2]]
  sign <- bf[[3]]
  writeRaster(type,paste0("changing_h_opt/classify/type/type_h",i,".tif"))
  writeRaster(subtype,paste0("changing_h_opt/classify/subtype/stype_h",i,".tif"))
  writeRaster(sign,paste0("changing_h_opt/classify/sign/sign_h",i,".tif"))
  rm(bf)
  #rm(out)
}

## Applying bfast01 to a raster to obtain breakpoint magnitude
bf01magn_raster<-function(matrix,raster,start,end,frequency,formula=response~trend+harmon,thrs=nlayers(raster)/2,h=0.15){
  r<-raster
  m<-matrix
  magn=c()
  for (x in 1:ncell(r)){
    if(sum(is.na(m[x,]))<thrs){  
      pix=c()
      for (i in 1:nlayers(r)){
        pix[i]<-as.vector(m[x,i])
      }
      ts<-ts(pix,start=start,end=end,frequency=frequency)
      b1<-bfast01(ts, formula = NULL, test = c("OLS-MOSUM", "BIC"), level = 0.05, aggregate = any,
                  trim = 0.2, bandwidth = h, functional = "max",
                  order = 0, lag = NULL, slag = NULL, na.action = na.omit, stl = "none")
      if(b1$breaks==1){
        b1.zoo <- as.zoo(b1) # data series
        ToB <- as.numeric(b1$breakpoints[[1]])  # time of break
        s1 <- b1$model[[2]]$coefficients[3] # slope segment 1
        s2 <- b1$model[[2]]$coefficients[4] # slope segment 2
        magn[x] <- as.numeric(b1.zoo$trend[ToB+1]) - as.numeric(b1.zoo$trend[ToB]) # magnitude of abrupt change
      }
    }else{
      magn[x]<-NA
    }
  }
  
  magn_raster<-raster(extent(r),nrow=nrow(r),ncol=ncol(r),crs=proj4string(r))
  values(magn_raster)<-magn
  
  return(magn_raster)
}

## Applying bfast01 to a raster to obtain breakpoint time
bf01time_raster<-function(matrix,raster,start,end,frequency,formula=response~trend+harmon,thrs=nlayers(raster)/2,h=0.15){
  r<-raster
  m<-matrix
  t=c()
  for (x in 1:ncell(r)){
    if(sum(is.na(m[x,]))<thrs){  
      pix=c()
      for (i in 1:nlayers(r)){
        pix[i]<-as.vector(m[x,i])
      }
      ts<-ts(pix,start=start,end=end,frequency=frequency)
      b1<-bfast01(ts, formula = NULL, test = c("OLS-MOSUM", "BIC"), level = 0.05, aggregate = any,
                  trim = 0.2, bandwidth = h, functional = "max",
                  order = 0, lag = NULL, slag = NULL, na.action = na.omit, stl = "none")
      
      if(b1$breaks==1){
        t[x]<-b1$breakpoints+1982 
      }
    }else{
      t[x]<-NA
    }
  }
  
  time_raster<-raster(extent(r),nrow=nrow(r),ncol=ncol(r),crs=proj4string(r))
  values(time_raster)<-t
  
  return(time_raster)
}


## Running bfast01 with different bandwidths
for(i in 2:33) {
  bf<-bf01c_raster(rue_crop,start=c(1982),end=c(2015),frequency=1,formula=response~trend,thrs=12,h=i/34)
  out<-stack(bf[[1]],bf[[2]],bf[[3]])
  writeRaster(out,paste0("changing_h_opt/classify/h",i,".tif"))
  rm(bf)
  rm(out)
}

## Running bfast01 with different bandwidths (to obtain time of break)
m<-as.matrix(rue_crop)
for(i in 2:33){
  tim_r<-bf01time_raster(m, rue_crop,start=c(1982),end=c(2015),frequency=1,formula=response~trend,thrs=12,h=i/34)
  writeRaster(tim_r,paste0("changing_h_opt/time/v2_time_h",i,".tif"))
  rm(tim_r)
}

## Running bfast01 with different bandwidths (to obtain magnitude of break)
m<-as.matrix(rue_crop)
for(i in 2:33){
  magn_r<-bf01magn_raster(m, rue_crop,start=c(1982),end=c(2015),frequency=1,formula=response~trend,thrs=12,h=i/34)
  writeRaster(magn_r,paste0("changing_h_opt/magn/v2_magn_h",i,".tif"))
  rm(tim_r)
}

#' Time Series Preprocessing for BFAST-Type Models
#' 
#' Time series preprocessing for subsequent regression modeling.  Based on a
#' (seasonal) time series, a data frame with the response, seasonal terms, a
#' trend term, (seasonal) autoregressive terms, and covariates is computed.
#' This can subsequently be employed in regression models.
#' 
#' To facilitate (linear) regression models of time series data, \code{bfastpp}
#' facilitates preprocessing and setting up regressor terms. It returns a
#' \code{data.frame} containing the first column of the \code{data} as the
#' \code{response} while further columns (if any) are used as covariates
#' \code{xreg}. Additionally, a linear trend, seasonal dummies, harmonic
#' seasonal terms, and (seasonal) autoregressive terms are provided.
#' 
#' Optionally, each column of \code{data} can be seasonally adjusted and/or
#' trend-adjusted via STL (season-trend decomposition via LOESS smoothing)
#' prior to preprocessing. The idea would be to capture season and/or trend
#' nonparametrically prior to regression modelling.
#' 
#' @param data A time series of class \code{\link[stats]{ts}}, or another
#' object that can be coerced to such. For seasonal components, a frequency
#' greater than 1 is required.
#' @param order numeric. Order of the harmonic term, defaulting to \code{3}.
#' @param lag numeric. Orders of the autoregressive term, by default omitted.
#' @param slag numeric. Orders of the seasonal autoregressive term, by default
#' omitted.
#' @param na.action function for handling \code{NA}s in the data (after all
#' other preprocessing).
#' @param stl character. Prior to all other preprocessing, STL (season-trend
#' decomposition via LOESS smoothing) can be employed for trend-adjustment
#' and/or season-adjustment.  The \code{"trend"} or \code{"seasonal"} component
#' or both from \code{\link[stats]{stl}} are removed from each column in
#' \code{data}. By default (\code{"none"}), no STL adjustment is used.
#' @param decomp "stlplus" or "stl": use the NA-tolerant decomposition package
#' or the reference package (which can make use of time series with 2-3
#' observations per year)
#' @param sbins numeric. Controls the number of seasonal dummies. If integer > 1,
#' sets the number of seasonal dummies to use per year.
#' If <= 1, treated as a multiplier to the number of observations per year, i.e.
#' `ndummies = nobs/year * sbins`.
#' @return If no formula is provided, \code{bfastpp} returns a
#' \code{"data.frame"} with the following variables (some of which may be
#' matrices).  \item{time}{numeric vector of time stamps,}
#' \item{response}{response vector (first column of \code{data}),}
#' \item{trend}{linear time trend (running from 1 to number of observations),}
#' \item{season}{factor indicating season period,} \item{harmon}{harmonic
#' seasonal terms (of specified \code{order}),} \item{lag}{autoregressive terms
#' (or orders \code{lag}, if any),} \item{slag}{seasonal autoregressive terms
#' (or orders \code{slag}, if any),} \item{xreg}{covariate regressor (all
#' columns of \code{data} except the first, if any).}
#' 
#' If a formula is given, \code{bfastpp} returns a \code{list} with components
#' \code{X}, \code{y}, and \code{t}, where \code{X} is the design matrix of the
#' model, \code{y} is the response vector, and \code{t} represents the time of
#' observations. \code{X} will only contain variables that occur in the
#' formula. Columns of \code{X} have names as decribed above.
#' @author Achim Zeileis
#' @seealso \code{\link[bfast]{bfastmonitor}}
#' @references \insertRef{janbfastmonitor}{bfast}
#' @keywords ts
#' @example examples/bfastpp.r
#' 
#' @export bfastpp
bfastpp<- function(data, order = 3,
                   lag = NULL, slag = NULL, na.action = na.omit,
                   stl = c("none", "trend", "seasonal", "both"),
                   decomp = c("stl", "stlplus"), sbins = 1)
{
  decomp <- match.arg(decomp)
  stl <- match.arg(stl)
  if (stl != "none" && decomp == "stlplus" && !requireNamespace("stlplus", quietly = TRUE))
    stop("Please install the stlplus package or set decomp = 'stl' or stl = 'none'.")
  
  ## double check what happens with 29-02 if that happens...
  ## we should keep it simple an remove the datum if that happens
  
  if(!is.ts(data)) data <- as.ts(data)
  
  ## STL pre-processing to try to adjust for trend or season
  stl <- match.arg(stl)
  if(stl != "none") {
    stl_adjust <- function(x) {
      x_stl <- if (decomp=="stlplus") {
        stlplus::stlplus(x, s.window = "periodic")$data
      } else {
        stats::stl(x, s.window = "periodic")$time.series
      }
      switch(stl,
             "trend" = x - x_stl[, "trend"],
             "seasonal" = x - x_stl[, "seasonal"],
             "both" = x - x_stl[, "trend"] - x_stl[, "seasonal"])
    }
    if(NCOL(data) > 1L) {
      for(i in 1:NCOL(data)) data[,i] <- stl_adjust(data[,i])
    } else {
      data <- stl_adjust(data)
    }
  }
  
  ## check for covariates
  if(NCOL(data) > 1L) {
    x <- coredata(data)[, -1L]
    y <- data[, 1L]
  } else {
    x <- NULL
    y <- data
  }
  
  ## data with trend and season factor
  seasonfreq <- if (sbins > 1) sbins else frequency(y)*sbins
  
  rval <- data.frame(
    time = as.numeric(time(y)),
    response = y,
    trend = 1:NROW(y),
    season = if (seasonfreq > 1) cut(cycle(y), seasonfreq, ordered_result = TRUE) else factor("no seasonality")
  )
  
  ## set up harmonic trend matrix as well
  freq <- frequency(y)
  order <- min(freq, order)
  harmon <- outer(2 * pi * as.vector(time(y)), 1:order)
  harmon <- cbind(apply(harmon, 2, cos), apply(harmon, 2, sin))
  colnames(harmon) <- if(order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == freq) harmon <- harmon[, -(2 * order)]
  rval$harmon <- harmon
  
  ## add lags
  nalag <- function(x, k) c(rep(NA, k), head(x, -k))
  if(!is.null(lag)) {
    rval$lag <- sapply(lag, function(k) nalag(as.vector(y), k))
    colnames(rval$lag) <- lag
  }
  if(!is.null(slag)) {
    rval$slag <- sapply(slag * freq, function(k) nalag(as.vector(y), k))
    colnames(rval$slag) <- slag
  }
  
  ## add regressors
  rval$xreg <- x
  
  ## omit missing values
  rval <- na.action(rval)
  
  ## return everything
  return(rval)
}

#' Checking for one major break in the time series
#' 
#' A function to select a suitable model for the data by choosing either a
#' model with 0 or with 1 breakpoint.
#' 
#' \code{bfast01} tries to select a suitable model for the data by choosing
#' either a model with 0 or with 1 breakpoint. It proceeds in the following
#' steps:
#' 
#' 1. The data is preprocessed with bfastpp using the arguments
#' `order`/`lag`/`slag`/`na.action`/`stl`/`sbins`.
#' 2. A linear model with the given formula is fitted. By default a suitable
#' formula is guessed based on the preprocessing parameters.
#' 3. The model with 1 breakpoint is estimated as well where the breakpoint is
#' chosen to minimize the segmented residual sum of squares.
#' 4. A sequence of tests for the null hypothesis of zero breaks is performed. Each
#' test results in a decision for FALSE (no breaks) or TRUE (structural
#' break(s)). The test decisions are then aggregated to a single decision (by
#' default using all() but any() or some other function could also be used).
#' 
#' Available methods for the object returned include standard methods for
#' linear models (coef, fitted, residuals, predict, AIC, BIC, logLik, deviance,
#' nobs, model.matrix, model.frame), standard methods for breakpoints
#' (breakpoints, breakdates), coercion to a zoo series with the decomposed
#' components (as.zoo), and a plot method which plots such a zoo series along
#' with the confidence interval (if the 1-break model is visualized). All
#' methods take a 'breaks' argument which can either be 0 or 1. By default the
#' value chosen based on the 'test' decisions is used.
#' 
#' Note that the different tests supported have power for different types of
#' alternatives. Some tests (such as supLM/supF or BIC) assess changes in all
#' coefficients of the model while residual-based tests (e.g., OLS-CUSUM or
#' OLS-MOSUM) assess changes in the conditional mean. See Zeileis (2005) for a
#' unifying view.
#' 
#' @param data A time series of class \code{\link[stats]{ts}}, or another
#' object that can be coerced to such. The time series is processed by
#' \code{\link[bfast]{bfastpp}}. A time series of class \code{\link[stats]{ts}}
#' can be prepared by a convenience function \code{\link[bfast]{bfastts}} in
#' case of daily, 10 or 16-daily time series.
#' @param formula formula for the regression model.  The default is
#' intelligently guessed based on the arguments order/lag/slag i.e.
#' \code{response ~ trend + harmon}, i.e., a linear trend and a harmonic season
#' component. Other specifications are possible using all terms set up by
#' \code{\link[bfast]{bfastpp}}, i.e., \code{season} (seasonal pattern with
#' dummy variables), \code{lag} (autoregressive terms), \code{slag} (seasonal
#' autoregressiv terms), or \code{xreg} (further covariates). See
#' \code{\link[bfast]{bfastpp}} for details.
#' @param test character specifying the type of test(s) performed. Can be one
#' or more of BIC, supLM, supF, OLS-MOSUM, ..., or any other test supported by
#' \code{\link[strucchangeRcpp]{sctest.formula}}
#' @param level numeric. Significance for the
#' \code{\link[strucchangeRcpp]{sctest.formula}} performed.
#' @param aggregate function that aggregates a logical vector to a single
#' value. This is used for aggregating the individual test decisions from
#' \code{test} to a single one.
#' @param trim numeric. The mimimal segment size passed to the \code{from}
#' argument of the \code{\link[strucchangeRcpp]{Fstats}} function.
#' @param bandwidth numeric scalar from interval (0,1), functional. The
#' \code{bandwidth} argument is passed to the \code{h} argument of the
#' \code{\link[strucchangeRcpp]{sctest.formula}}.
#' @param functional arguments passed on to
#' \code{\link[strucchangeRcpp]{sctest.formula}}
#' @param order numeric. Order of the harmonic term, defaulting to \code{3}.
#' @param lag numeric. Order of the autoregressive term, by default omitted.
#' @param slag numeric. Order of the seasonal autoregressive term, by default
#' omitted.
#' @param na.action arguments passed on to \code{\link[bfast]{bfastpp}}
#' @param reg whether to use OLS regression [lm()] or
#' robust regression [MASS::rlm()].
#' @param stl argument passed on to \code{\link[bfast]{bfastpp}}
#' @param sbins argument passed on to \code{\link[bfast]{bfastpp}}
#' @return \code{bfast01} returns a list of class \code{"bfast01"} with the
#' following elements: \item{call}{the original function call.} \item{data}{the
#' data preprocessed by \code{"bfastpp"}.} \item{formula}{the model formulae.}
#' \item{breaks}{the number of breaks chosen based on the \code{test} decision
#' (either 0 or 1).} \item{test}{the individual test decisions.}
#' \item{breakpoints}{the optimal breakpoint for the model with 1 break.}
#' \item{model}{A list of two 'lm' objects with no and one breaks,
#' respectively.}
#' @author Achim Zeileis, Jan Verbesselt
#' @seealso \code{\link[bfast]{bfastmonitor}},
#' \code{\link[strucchangeRcpp]{breakpoints}}
#' @references \insertRef{rogierbfast01}{bfast}
#'
#' \insertRef{achimstrucchange}{bfast}
#' @keywords ts
#' @examples
#' 
#' library(zoo)
#' ## define a regular time series
#' ndvi <- as.ts(zoo(som$NDVI.a, som$Time))
#' 
#' ## fit variations
#' bf1 <- bfast01(ndvi)
#' bf2 <- bfast01(ndvi, test = c("BIC", "OLS-MOSUM", "supLM"), aggregate = any)
#' bf3 <- bfast01(ndvi, test = c("OLS-MOSUM", "supLM"), aggregate = any, bandwidth = 0.11) 
#' 
#' ## inspect test decisions
#' bf1$test
#' bf1$breaks
#' bf2$test
#' bf2$breaks
#' bf3$test
#' bf3$breaks
#' 
#' ## look at coefficients
#' coef(bf1)
#' coef(bf1, breaks = 0)
#' coef(bf1, breaks = 1) 
#' 
#' ## zoo series with all components
#' plot(as.zoo(ndvi))
#' plot(as.zoo(bf1, breaks = 1))
#' plot(as.zoo(bf2))
#' plot(as.zoo(bf3))
#' 
#' ## leveraged by plot method
#' plot(bf1, regular = TRUE)
#' plot(bf2)
#' plot(bf2, plot.type = "multiple",
#'      which = c("response", "trend", "season"), screens = c(1, 1, 2))
#' plot(bf3)
#' 
#' 
#' @export bfast01
bfast01 <- function(data, formula = NULL,
                    test = "OLS-MOSUM", level = 0.05, aggregate = all,
                    trim = NULL, bandwidth = 0.15, functional = "max",
                    order = 3, lag = NULL, slag = NULL, na.action = na.omit, 
                    reg = c("lm", "rlm"), stl = "none", sbins = 1)
{
  reg <- match.arg(reg)
  if (reg == "rlm"){
    if (!requireNamespace("MASS"))
      stop("reg = 'rlm' requires the MASS package")
    reg <- MASS::rlm
  }
  ## data preprocessing
  stl <- match.arg(stl, c("none", "trend", "seasonal", "both"))
  if (!inherits(data, "data.frame")) 
    data <- bfastpp(data, order = order, lag = lag, slag = slag, 
                    na.action = na.action, stl = stl, sbins = sbins)
  if (is.null(formula)) {
    formula <- c(trend = !(stl %in% c("trend", "both")), 
                 harmon = order > 0 & !(stl %in% c("seasonal", "both")), 
                 lag = !is.null(lag), slag = !is.null(slag))
    formula <- as.formula(paste("response ~", paste(names(formula)[formula], 
                                                    collapse = " + ")))
  }
  ## fit 1-segment model
  model1 <- do.call(reg, list(formula = formula, data = data))
  ## determine optimal single breakpoint
  if(is.null(trim)) trim <- 5 * length(coef(model1))
  fs <- Fstats(formula, data = data, from = trim)
  bp <- breakpoints(fs)
  ## fit 2-segment model
  data$segment <- breakfactor(bp)
  levels(data$segment) <- c("1", "2")
  formula2 <- update(update(formula, . ~ segment/(.)), . ~ . - 1)
  model2   <- do.call(reg, list(formula = formula2, data = data) )
  ## compute BIC values
  bic <- c(BIC(model1), BIC(model2) + log(nrow(data)))
  ## perform tests
  improvement01 <- function(test) {
    trim01 <- if(trim > 1) trim/nrow(data) else trim
    if (test == "BIC") return(bic[2] < bic[1])
    if (test %in% c("supF", "aveF", "expF")) 
      return( sctest(fs, type = test)$p.value < level)
    if (test == "supLM") 
      return( sctest(gefp(formula, data = data), 
                     functional = supLM(trim01))$p.value < level )
    sctest(formula, data = data, type = test, h = bandwidth, 
           functional = functional)$p.value < level
  }
  test <- structure(sapply(test, improvement01), names = test)
  rval <- list(call = match.call(), data = data, formula = formula, 
               breaks = as.numeric(aggregate(test)), breakpoints = bp$breakpoints,
               # DM: Could return 0 if the breakpoint is insignificant using breakpoints=ifelse(as.numeric(aggregate(test))==0,0,bp$breakpoints)
               test = test, model = list(model1, model2))
  class(rval) <- "bfast01"
  rval$confint <- .confint01(rval, level = 1 - level)
  return(rval)
}

#' @method breakpoints bfast01
#' @export
breakpoints.bfast01 <- function(obj, breaks = NULL, ...) {
  if(is.null(breaks)) breaks <- obj$breaks
  n <- nrow(obj$data)
  rval <- list(
    breakpoints = if(breaks > 0) obj$breakpoints else NA,
    RSS = if(breaks > 0) deviance(obj$model[[1]]) else deviance(obj$model[[2]]),
    nobs = n,
    nreg = length(coef(obj$model[[1]])),
    call = match.call(),
    datatsp = c(1/n, 1, n)
  )
  class(rval) <- "breakpoints"
  return(rval)
}

#' @method breakdates bfast01
#' @export
breakdates.bfast01 <- function(obj, format.times = NULL, breaks = NULL, ...) {
  if(is.null(breaks)) breaks <- obj$breaks
  if(breaks > 0) obj$data$time[obj$breakpoints] else NA
}

#' @method logLik bfast01
#' @export
logLik.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  rval <- logLik(object$model[[breaks + 1]])
  attr(rval, "df") <- attr(rval, "df") + breaks
  rval
}

#' @method deviance bfast01
#' @export
deviance.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  deviance(object$model[[breaks + 1]])
}

#' @method model.frame bfast01
#' @export
model.frame.bfast01 <- function(formula, breaks = NULL, ...) model.frame(formula$model[[1]])

#' @method model.matrix bfast01
#' @export
model.matrix.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  model.matrix(object$model[[breaks + 1]])
}

#' @method nobs bfast01
#' @export
nobs.bfast01 <- function(object, breaks = NULL, ...) nrow(object$data)

#' @method AIC bfast01
#' @export
AIC.bfast01 <- function(object, breaks = NULL, ...) AIC(logLik(object, breaks = breaks), ...)
#' @method BIC bfast01
#' @export
BIC.bfast01 <- function(object, breaks = NULL, ...) BIC(logLik(object, breaks = breaks), ...)

#' @method coef bfast01
#' @export
coef.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  cf0 <- coef(object$model[[1]])
  if(breaks < 1) return(cf0)
  cf <- matrix(coef(object$model[[2]]), nrow = 2)
  colnames(cf) <- names(cf0)
  bd <- object$data$time[c(1, object$breakpoints, object$breakpoints + 1, nrow(object$data))]
  bd <- format(round(bd, digits = 3))
  rownames(cf) <- paste(bd[c(1, 3)], bd[c(2, 4)], sep = "--")
  cf
}

#' @method fitted bfast01
#' @export
fitted.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  fitted(object$model[[breaks + 1]])
}

#' @method residuals bfast01
#' @export
residuals.bfast01 <- function(object, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  residuals(object$model[[breaks + 1]])
}

#' @method predict bfast01
#' @export
predict.bfast01 <- function(object, newdata, breaks = NULL, ...) {
  breaks <- .breaks01(object, breaks)
  predict(object$model[[breaks + 1]], newdata, ...)
}

#' @method as.zoo bfast01
#' @export
as.zoo.bfast01 <- function(x, breaks = NULL, ...) {
  breaks <- .breaks01(x, breaks)
  
  ## fitted values
  d <- x$data
  fit <- predict(x, newdata = d, breaks = breaks)
  
  ## residuals
  res <- x$data$response - fit
  
  ## eliminate seasonal effects
  if(!is.null(d$harmon)) d$harmon <- d$harmon * 0
  if(!is.null(d$season)) d$season <- levels(d$season)[1]
  season <- fit - predict(x, newdata = d, breaks = breaks)
  
  ## eliminate (auto)regressive effects
  for(i in c("lag", "slag", "xreg")) if(!is.null(d[[i]])) d[[i]] <- d[[i]] * 0
  reg <- fit - season - predict(x, newdata = d, breaks = breaks)
  
  ## compute fit = trend + season + reg
  trend <- fit - season - reg
  
  ## include mean in trend instead of reg
  m <- if(breaks > 0) tapply(reg, x$data$segment, mean)[x$data$segment] else mean(reg)
  trend <- trend + m
  reg <- reg - m
  
  rval <- cbind(x$data$response, fit, trend, season, reg, res)
  colnames(rval) <- c("response", "fitted", "trend", "season", "reg", "residuals")
  zoo(rval, x$data$time)
}

#' @method plot bfast01
#' @export
plot.bfast01 <- function(x, breaks = NULL, which = c("response", "fitted", "trend"),
                         plot.type = "single", panel = NULL, screens = NULL,
                         col = NULL, lwd = NULL,
                         main = "", xlab = "Time", ylab = NULL, ci = NULL, regular = TRUE, ...)
{
  ## set up zoo series and select series to be plotted
  breaks <- .breaks01(x, breaks)
  z <- as.zoo(x, breaks = breaks)
  which <- sapply(which, function(x) match.arg(x, colnames(z)))
  z <- z[, which]
  
  ## try making intelligent guesses about default col/lwd
  plot.type <- match.arg(plot.type, c("single", "multiple"))
  if(is.null(col)) {
    col0 <- c("gray", "black", "blue", "red", "green", "black")
    if(plot.type == "single") {
      col <- col0
      names(col) <- c("response", "fitted", "trend", "season", "reg", "residuals")
      col <- col[which]
    } else {
      col <- if(is.null(screens)) 1 else unlist(lapply(unique(screens),
                                                       function(i) if((n <- sum(screens == i)) == 1) "black" else rep(col0, length.out = n)))
    }
  }
  if(is.null(lwd)) {
    if(plot.type == "single") {
      lwd <- c(2, 1, 2, 1, 1, 1)
      names(lwd) <- c("response", "fitted", "trend", "season", "reg", "residuals")
      lwd <- lwd[which]
    } else {
      lwd <- 1
    }
  }
  
  ## default y-axis labels
  if(is.null(ylab)) {
    ylab <- which
    for(i in seq_along(ylab)) substr(ylab[i], 1, 1) <- toupper(substr(ylab[i], 1, 1))
    if(plot.type == "single") {
      ylab <- paste(ylab, collapse = " / ")
    } else {
      if(!is.null(screens)) ylab <- sapply(unique(screens), function(i) paste(ylab[screens == i], collapse = " / "))
    }
  }
  
  ## set up panel function with confidence intervals
  if(is.null(panel)) panel <- .make_confint_lines01(x, breaks = breaks, ci = ci)
  
  if(regular) z <- as.zoo(as.ts(z))
  
  plot(z, plot.type = plot.type, panel = panel, screens = screens,
       col = col, lwd = lwd, main = main, xlab = xlab, ylab = ylab, ...)
}

.breaks01 <- function(object, breaks) {
  if(is.null(breaks)) breaks <- object$breaks
  breaks <- breaks[1]
  if(!breaks %in% 0:1) stop("breaks can only be 0 or 1")
  breaks
}

.make_confint_lines01 <- function(object, breaks = NULL, col = 1, lty = 2, lwd = 1, ci = list(), ...)
{
  breaks <- .breaks01(object, breaks)
  if(breaks < 1) return(lines)
  
  function(x, y, ...) {
    lines(x, y, ...)
    abline(v = breakdates(object, breaks = breaks), lty = lty, col = col, lwd = lwd)
    if(!identical(ci, FALSE)) {
      if(!is.list(ci)) ci <- list()
      if(is.null(ci$col)) ci$col <- 2
      if(is.null(ci$angle)) ci$angle <- 90
      if(is.null(ci$length)) ci$length <- 0.05
      if(is.null(ci$code)) ci$code <- 3
      if(is.null(ci$at)) {
        at <- par("usr")[3:4]
        at <- diff(at)/1.08 * 0.02 + at[1]
      } else {
        at <- ci$at
      }
      ci$at <- NULL
      do.call("arrows", c(list(
        x0 = object$data$time[object$confint[1]],
        y0 = at,
        x1 = object$data$time[object$confint[3]],
        y1 = at),
        ci))
    }
  }
}

.confint01 <- function(object, level = 0.95, het.reg = TRUE, het.err = TRUE)
{
  ## data and arguments
  X <- model.matrix(object$model[[1]])
  y <- model.response(model.frame(object$model[[1]]))
  n <- nrow(object$data)
  a2 <- (1 - level)/2
  bp <- c(0, object$breakpoints, n)
  
  ## auxiliary functions
  myfun <- function(x, level = 0.975, xi = 1, phi1 = 1, phi2 = 1)
    (pargmaxV(x, xi = xi, phi1 = phi1, phi2 = phi2) - level)
  myprod <- function(delta, mat) as.vector(crossprod(delta, mat) %*% delta)
  
  ## overall fits
  res <- residuals(object$model[[2]])
  beta <- coef(object, breaks = 1)
  sigma1 <- sigma2 <- sum(res^2)/n
  Q1 <- Q2 <- crossprod(X)/n
  Omega1 <- Omega2 <- sigma1 * Q1
  
  ## subsample fits
  X1 <- X[(bp[1]+1):bp[2],,drop = FALSE]
  X2 <- X[(bp[2]+1):bp[3],,drop = FALSE]
  y1 <- y[(bp[1]+1):bp[2]]
  y2 <- y[(bp[2]+1):bp[3]]
  beta1 <- beta[1,]
  beta2 <- beta[2,]
  if(het.reg) {
    Q1 <- crossprod(X1)/nrow(X1)
    Q2 <- crossprod(X2)/nrow(X2)
  }
  if(het.err) {
    sigma1 <- sum(res[(bp[1]+1):(bp[2])]^2)/nrow(X1)
    sigma2 <- sum(res[(bp[2]+1):(bp[3])]^2)/nrow(X2)
    Omega1 <- sigma1 * Q1
    Omega2 <- sigma2 * Q2
  }
  delta <- beta2 - beta1
  
  Oprod1 <- myprod(delta, Omega1)
  Oprod2 <- myprod(delta, Omega2)
  Qprod1 <- myprod(delta, Q1)
  Qprod2 <- myprod(delta, Q2)
  
  xi <- if(het.reg) Qprod2/Qprod1 else 1
  phi1 <- sqrt(Oprod1/Qprod1)
  phi2 <- sqrt(Oprod2/Qprod2)
  
  p0 <- pargmaxV(0, phi1 = phi1, phi2 = phi2, xi = xi)
  if(is.nan(p0) || p0 < a2 || p0 > (1-a2)) {
    warning(paste("Confidence interval cannot be computed: P(argmax V <= 0) =", round(p0, digits = 4)))
    upper <- NA
    lower <- NA
  } else {
    ub <- lb <- 0
    while(pargmaxV(ub, phi1 = phi1, phi2 = phi2, xi = xi) < (1 - a2)) ub <- ub + 1000
    while(pargmaxV(lb, phi1 = phi1, phi2 = phi2, xi = xi) > a2) lb <- lb - 1000
    
    upper <- uniroot(myfun, c(0, ub), level = (1-a2), xi = xi, phi1 = phi1, phi2 = phi2)$root
    lower <- uniroot(myfun, c(lb, 0), level = a2, xi = xi, phi1 = phi1, phi2 = phi2)$root
    
    upper <- upper * phi1^2 / Qprod1
    lower <- lower * phi1^2 / Qprod1
  }
  
  bp <- c(bp[2] - ceiling(upper), bp[2], bp[2] - floor(lower))
  a2 <- round(a2 * 100, digits = 1)
  names(bp) <- c(paste(a2, "%"), "breakpoints", paste(100 - a2, "%"))
  bp
}

