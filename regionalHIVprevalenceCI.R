### title: "BIOST 555 Project"
### author: "Yao He"
### date: "03/18/2020"

# Load packages
library(knitr); library(rgdal); library(maptools); library(sp); library(spdep); library(SpatialEpi); 
library(RColorBrewer); library(ggplot2); library(maps); library(broom); library(INLA); library(SUMMER); library(survey);
library(rgeos); library(gpclib); library(haven)

## read in data
hiv <- read_dta("CIAR61FL.DTA") # read in hiv data; as tbl_df which will cause problem w/i fitGeneric fx
hiv <- as.data.frame(hiv) # convert the tbl_df to data.frame
hiv$hhid <- paste(hiv$hivclust,"-",hiv$hivnumb) # create household unique id
hiv2 <- hiv[-which(hiv$hiv03==7),] # drop the undetermined test results; there were 4 people
ci_map2 <- readOGR(dsn = ".", layer = "CIGE61FL") # read in GPS data
ci_regions <- readOGR(dsn = ".", layer = "sdr_subnational_boundaries") # read in regional boundaries

# add urban/rural, region names & codes to HIV dataset
hiv2$urban <- ci_map2@data$URBAN_RURA[match(hiv2$hivclust,ci_map2@data$DHSCLUST)]
hiv2$regna <- ci_map2@data$DHSREGNA[match(hiv2$hivclust,ci_map2@data$DHSCLUST)]
hiv2$regna <- as.character(hiv2$regna)
hiv2$regco <- ci_map2@data$DHSREGCO[match(hiv2$hivclust,ci_map2@data$DHSCLUST)]

# add country code to HIV dataset
hiv2$country <- 1

# create strata variable based on survey design
hiv2$strata <- with(hiv2, ifelse(urban=="U", regco, regco+11))
hiv2$strata <- as.integer(hiv2$strata)

# count # of households
length(unique(hiv$hhid)) # in initial dataset
length(unique(hiv2$hhid)) # in final dataset

## table 1 characteristics
# test result
table(hiv2$hiv03)
# urban vs rural
tab.ur <- with(hiv2, table(urban,hiv03)) 
tab.ur
round(prop.table(tab.ur,2)*100,1)
tab.ur1 <- with(hiv2, table(urban)) 
tab.ur1
round(prop.table(tab.ur1)*100,1)
# by region
tab.reg <- with(hiv2, table(regna,hiv03))
tab.reg
round(prop.table(tab.reg,2)*100,1)
tab.reg1 <- with(hiv2, table(regna)) 
tab.reg1
round(prop.table(tab.reg1)*100,1)

# adjacency matrix
nb.r2 <- poly2nb(ci_regions, queen = F, row.names = ci_regions$DHSREGFR)
mat2 <- nb2mat(nb.r2, style = "B", zero.policy = TRUE)
colnames(mat2) <- rownames(mat2)
mat2 <- as.matrix(mat2[1:dim(mat2)[1], 1:dim(mat2)[1]])

# Weighted (direct) estimates for each region
design <- svydesign(ids = ~1, weights = ~hiv05, 
                    strata = ~strata, data = hiv2)
tab.w.0 <- svyby(~hiv03, ~country, design, svymean)
colnames(tab.w.0) <- c("region","weighted","se")
tab.w.1 <- svyby(~hiv03, ~regna, design, svymean)
colnames(tab.w.1) <- c("region","weighted","se")
tab.w <- rbind(tab.w.0, tab.w.1)
tab.w$lower <- tab.w$weighted - 1.96*tab.w$se
tab.w$higher <- tab.w$weighted + 1.96*tab.w$se
tab.w[,-1]

# Smoothed weighted direct estimates (again take the posterior medians).**
svysmoothed <- fitGeneric(data = hiv2, geo = ci_regions, Amat = mat2, 
                          responseType = "binary",responseVar = "hiv03", 
                          strataVar = "strata", weightVar = "hiv05",
                          regionVar = "regna", clusterVar = "~hivclust+hivnumb", CI = 0.95)
                          
## table of smoothed posterior medians
tab.w.sm1 <- with(svysmoothed, cbind(smooth$median.original, 
                                       sqrt(smooth$variance.original)))
colnames(tab.w.sm1) <- c("w.sm.median","sd")
rownames(tab.w.sm1) <- rownames(mat2)
tab.w.sm1

## table of smoothed posterior means
tab.w.sm2 <- with(svysmoothed, cbind(smooth$mean.original,
                                     smooth$lower.original,
                                     smooth$upper.original))
colnames(tab.w.sm2) <- c("w.sm.mean","lower","upper")
rownames(tab.w.sm2) <- rownames(mat2)
tab.w.sm2

## plots
toplot <- svysmoothed$smooth
toplot$HTest <- svysmoothed$HT$HT.est.original
lo <- svysmoothed$HT$HT.est.original-1.96*sqrt(svysmoothed$HT$HT.variance.original)
hi <- svysmoothed$HT$HT.est.original+1.96*sqrt(svysmoothed$HT$HT.variance.original)
toplot$HTlower <- lo
toplot$HTupper <- hi

mapPlot(data=toplot, geo=ci_regions,
        variables=c("mean.original","median.original"),
        labels=c("Smoothed posterior mean","Smoothed posterior median"),
        by.data="region", by.geo="DHSREGFR")

mapPlot(data=toplot, geo=ci_regions,
        variables=c("HTest","HTlower","HTupper",
                    "mean.original","lower.original","upper.original"),
        labels=c("Direct estimates","Direct lower","Direct upper",
                 "Smoothed posterior mean","Smoothed lower","Smoothed upper"),
        by.data="region", by.geo="DHSREGFR")

### A lot more variability in the direct weighted estimates than the smoothed weighted estimates.

# rank weighted estimates
wsrank.w <- tab.w
wsrank.w$rank <- rank(-wsrank.w$weighted)
head(wsrank.w[order(wsrank.w$rank),], 12)

# rank smoothed estimates
wsrank.s <- as.data.frame(tab.w.sm2)
wsrank.s$rank <- rank(-wsrank.s$w.sm.mean)
head(wsrank.s[order(wsrank.s$rank),], 12)

# compile point estimates for comprison plots
est <- data.frame(weighted = svysmoothed$HT$HT.est.original,
                  weightedsmooth = svysmoothed$smooth$mean.original)
# comparison plots for estimate
l1 <- range(est)
plot(weightedsmooth ~ weighted, data = est,
     pch = 19, col = "skyblue3", xlim = l1, ylim = l1,
     ylab = "Smoothed weighted posterior mean", xlab = "Weighted direct estimates")
abline(0, 1, col = "gray15")

# compare standard deviations
toplot$weighted.sd <- sqrt(svysmoothed$HT$HT.variance.original)
toplot$posterior.sd <- sqrt(toplot$variance.original)
sd <- data.frame(weighted.sd = toplot$weighted.sd,
                 weightedsmooth.sd = toplot$posterior.sd)
l2 <- range(sd)
plot(weightedsmooth.sd ~ weighted.sd, data = sd,
     pch = 19, col = "deeppink", xlim = l2, ylim = l2,
     ylab = "Smoothed weighted posterior SD", xlab = "Weighted direct SD")
abline(0, 1, col = "gray15")
