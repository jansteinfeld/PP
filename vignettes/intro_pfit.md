<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Getting started with Personfit in PP}
-->
# Getting started with Person-Fit functions

A brief introduction in the currently implemented person fit functions will be added soon. Currently the **LZ, LZ*** and also the **Infit-Outfit-Statistics** are implemented. We also added the Infit-Outfit-Functions for the Partial-Credit Model. Meanwhile we are working on plots for a better understanding of the person misfit as well as on inference statistic methods.

Now a simple example.
First we will simulate some data for our hands on example:

```r
library(PP)

set.seed(1337)

# simulate some intercepts
diffpar <- seq(-3,3,length=15)
# simulate some slope parameters
sl     <- round(runif(15,0.5,1.5),2)
la     <- round(runif(15,0,0.25),2)
ua     <- round(runif(15,0.8,1),2)

# simulate response matrix
awm <- matrix(sample(0:1,100*15,replace=TRUE),ncol=15)
```

1. We will start with a simple 1PL-Model. First we have to estimate the person parameters. Here we have to choose an estimation method. It is importend, that you can only choose **mle, wle or map** for the LZ and LZ* Index. For the Infit-Outfit statistic we support the **mle and wle** estimates.


```r
# MLE
res1plmle <- PP_4pl(respm = awm,thres = diffpar,type = "mle")
```

```
## Estimating:  1pl model ... 
## type = mle 
## Estimation finished!
```

```r
# WLE
res1plwle <- PP_4pl(respm = awm,thres = diffpar,type = "wle")
```

```
## Estimating:  1pl model ... 
## type = wle 
## Estimation finished!
```

```r
# MAP estimation
res1plmap <- PP_4pl(respm = awm,thres = diffpar,type = "map")
```

```
## Warning in PP_4pl(respm = awm, thres = diffpar, type = "map"): all mu's are set to 0!
```

```
## Warning in PP_4pl(respm = awm, thres = diffpar, type = "map"): all sigma2's are set to 1!
```

```
## Estimating:  1pl model ... 
## type = map 
## Estimation finished!
```

We also support the 2PL, 3PL and 4PL Model:

```r
# ------------------------------------------------------------------------
## 2PL model ##### 
# ------------------------------------------------------------------------
# MLE
res2plmle <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,type = "mle")
```

```
## Estimating:  2pl model ... 
## type = mle 
## Estimation finished!
```

```r
# WLE
res2plwle <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,type = "wle")
```

```
## Estimating:  2pl model ... 
## type = wle 
## Estimation finished!
```

```r
# ------------------------------------------------------------------------
## 3PL model ##### 
# ------------------------------------------------------------------------
# MLE
res3plmle <- PP_4pl(respm = awm,thres = diffpar,
                    slopes = sl,lowerA = la,type = "mle")
```

```
## Estimating:  3pl model ... 
## type = mle 
## Estimation finished!
```

```r
# WLE
res3plwle <- PP_4pl(respm = awm,thres = diffpar,
                    slopes = sl,lowerA = la,type = "wle")
```

```
## Estimating:  3pl model ... 
## type = wle 
## Estimation finished!
```

```r
# ------------------------------------------------------------------------
## 4PL model ##### 
# ------------------------------------------------------------------------
# MLE
res4plmle <- PP_4pl(respm = awm,thres = diffpar,
                    slopes = sl,lowerA = la,upperA=ua,type = "mle")
```

```
## Estimating:  4pl model ... 
## type = mle 
## Estimation finished!
```

```r
# WLE
res4plwle <- PP_4pl(respm = awm,thres = diffpar,
                    slopes = sl,lowerA = la,upperA=ua,type = "wle")
```

```
## Estimating:  4pl model ... 
## type = wle 
## Estimation finished!
```

2. After the estimation of the person parameter we are able to calculate the personfits. At this point you are able to calculate only one kind of personfit as well as all simultaneously (as shown next).


```r
# ------------------------------------------------------------------------
## 1PL model ##### 
# ------------------------------------------------------------------------
## LZ*-Index ##### 
pfit1pl_lz <- Pfit(respm=awm,pp=res1plwle,fitindices="lzstar")
## LZ*-Index combined with Infit-Outfit ##### 
pfit1pl_li <- Pfit(respm=awm,pp=res1plwle,fitindices=c("lzstar","infitoutfit"))
# ------------------------------------------------------------------------
## 2PL model ##### 
# ------------------------------------------------------------------------
## LZ*-Index ##### 
pfit2pl_lz <- Pfit(respm=awm,pp=res2plwle,fitindices="lzstar")
## LZ*-Index combined with Infit-Outfit ##### 
pfit2pl_li <- Pfit(respm=awm,pp=res2plwle,fitindices=c("lzstar","infitoutfit"))
# ------------------------------------------------------------------------
## 3PL model ##### 
# ------------------------------------------------------------------------
## LZ*-Index ##### 
pfit3pl_lz <- Pfit(respm=awm,pp=res3plwle,fitindices="lzstar")
## LZ*-Index combined with Infit-Outfit ##### 
pfit3pl_li <- Pfit(respm=awm,pp=res3plwle,fitindices=c("lzstar","infitoutfit"))
# ------------------------------------------------------------------------
## 4PL model ##### 
# ------------------------------------------------------------------------
## LZ*-Index ##### 
pfit4pl_lz <- Pfit(respm=awm,pp=res4plwle,fitindices="lzstar")
## LZ*-Index combined with Infit-Outfit ##### 
pfit4pl_li <- Pfit(respm=awm,pp=res4plwle,fitindices=c("lzstar","infitoutfit"))
```

3. We can also use different person parameter estimates

```r
# ------------------------------------------------------------------------
## 1PL model ##### 
# ------------------------------------------------------------------------
## LZ*-Index ##### 
## mle ####
pfit1pl_mle_l <- Pfit(respm=awm,pp=res1plmle,fitindices="lzstar")
## wle ####
pfit1pl_wle_l <- Pfit(respm=awm,pp=res1plwle,fitindices="lzstar")
## map ####
pfit1pl_map_l <- Pfit(respm=awm,pp=res1plmap,fitindices="lzstar")
```