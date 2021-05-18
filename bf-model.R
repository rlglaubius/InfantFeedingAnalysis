## Robert Glaubius
## Avenir Health
## 2021-05-18
##
## bf-model.R uses data from household surveys with HIV testing done in
## sub-Saharan Africa to estimate breastfeeding duration among mothers
## stratified by HIV status

## Load packages
library(data.table)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(rdhs)
library(readxl)
library(rstan)

source("bf-support.R")

## +=+ Function definitions +==================================================+

## Breastfeeding model
## months - months since last birth
## med  - month at which 50% of mothers have stopped breastfeeding among those who breastfed initially
## shp  - log-logistic shape parameter
## init - proportion of mothers who breastfed initially
bf.model = function(months, med, shp, init) {
  return(init * (1.0 - 1.0 / (1.0 + (months / med)^-shp)))
}

## Return the median duration of breastfeeding (the months since birth at which
## half of women are no longer breastfeeding)
bf.median = function(med, shp, init) {
  return(med * (1/(1 - 0.5/init) - 1)^(-1/shp))
}

## Fit the model to HIV-negative women in a single survey
## data.flat - data frame of breastfeeding estimates
## Return the stanfit object returned by stan() (ssee ?stan for details)
bf.fit = function(data.flat) {
  ## extract a data frame that has one row per survey in the order they occur in data.flat
  abstract = data.flat[!duplicated(data.flat[,c('SubregionName', 'CountryName', 'SurveyYear', 'SurveyId')]),]

  months = c(1.5, seq(4,36,2)) # midpoints of age intervals used for PHIA datasets [0,2] (2,4] (4,6] (6,8] ... (34,36]
  # months = seq(1, 35, 2) # midpoints of Spectrum's age intervals [0,2) [2,4) [4,6) [6,8) ... [34,36)
  
  n.months = length(months)
  n.survey = nrow(abstract)
  n.country = length(unique(abstract$CountryName))
  n.region = length(unique(abstract$SubregionName))
  
  survey.dat = list(
    n_surveys = n.survey,
    n_countries = n.country,
    n_regions = n.region,
    n_months = n.months,
    month = months,
    pmtct = abstract$PMTCT_coverage,
    country = as.array(as.integer(factor(abstract$CountryName))),
    region = as.array(as.integer(factor(abstract$SubregionName))),
    survey_year = as.numeric(abstract$SurveyYear),
    ref_year = 2010,
    X = array(data.flat$not.bf, dim=c(n.months, n.survey, 2)),
    Y = array(data.flat$yes.bf, dim=c(n.months, n.survey, 2)))

  param.init = function(chain_id=1) {
    list(
      base_bgn_hiv = as.array(rep(0, n.region)),
      base_med_hiv = as.array(rep(0, n.region)),

      base_bgn_arv = as.array(rep(0, n.region)),
      base_med_arv = as.array(rep(0, n.region)),

      base_slope_bgn = as.array(rep(0, n.region)),
      base_slope_med = as.array(rep(0, n.region)),
      base_slope_shp = as.array(rep(0, n.region)),

      ceff_bgn = as.array(rep(0.8,     n.country)),
      ceff_med = as.array(rep(log(22), n.country)),
      ceff_shp = as.array(rep(log(4),  n.country))
    )
  }
  
  fit = stan(file = "bf-model.stan",
             data = survey.dat,
             init = param.init,
             iter = 2000,
             chains = 4)
  return(fit)
}

plot.model.fit = function(stan.fit, data.flat, plot.file) {
  abstract = data.flat[!duplicated(data.flat[,c('SubregionName', 'CountryName', 'SurveyYear', 'SurveyId')]),]

  months = c(1.5, seq(4,36,2)) # midpoints of age intervals used for PHIA datasets [0,2] (2,4] (4,6] (6,8] ... (34,36]
  # months = seq(1, 35, 2) # midpoints of Spectrum's age intervals [0,2) [2,4) [4,6) [6,8) ... [34,36)
  surveyIds = abstract$SurveyId
  n.survey = length(surveyIds)

  draws = as.data.frame(stan.fit)
  fit.summary = summary(stan.fit)$summary
  post.max = which.max(draws$lp__)
  post.ind = array(grep('propbf', colnames(draws)), dim=c(length(months), n.survey, 2))
  summ.ind = array(grep('propbf', rownames(fit.summary)), dim=c(2, n.survey, length(months)))
  
  gset = runif(nrow(draws), 0.6, 0.8)
  cset = rgb(gset, gset, gset)
  
  for (s in 1:n.survey) {
    sid = which(abstract$SurveyId == surveyIds[s])

    plot.file.s = sprintf("Figures/%s-%s.tiff", gsub(".tiff", "", plot.file), surveyIds[s])
    cat(sprintf("Plotting %s (%s %s %s)\n", plot.file.s, abstract$CountryName[sid], abstract$SurveyYearLabel[sid], abstract$SurveyType[sid]))

    data.neg = data.flat[which(data.flat$SurveyId==surveyIds[s] & data.flat$hiv=='negative'),]
    data.pos = data.flat[which(data.flat$SurveyId==surveyIds[s] & data.flat$hiv=='positive'),]

    tiff(plot.file.s, w=3*3.42, h=2.44, units="in", pointsize=8, compression="lzw", res=300)
    layout(matrix(1:3, nrow=1, byrow=TRUE))
    
    ## 1. Fit to HIV- women, with data
    par(las=1, mar=c(2.5,3.0,1.5,1.0), cex=1)
    plot(c(0,36), c(0,1), cex=0, ann=FALSE, axes=FALSE, yaxs='i')
    # for (k in 1:nrow(draws)) {lines(months, draws[k,post.ind[,s,1]], col=cset[k])}
    polygon(c(months, rev(months)), c(fit.summary[summ.ind[1,s,],4], rev(fit.summary[summ.ind[1,s,],8])), col='#0066ff30', lty=0)
    lines(months, draws[post.max,post.ind[,s,1]], col='#0066ff')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=seq(0,36,6), labels=seq(0,36,6))
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    arrows(months[data.neg$age], data.neg$lower, months[data.neg$age], data.neg$upper, code=3, angle=90, length=0.0125, col='#000000', xpd=1)
    points(months[data.neg$age], data.neg$value, pch=16, col='#000000', xpd=1)
    legend('bottomleft', inset=c(0.025,0), box.lty=0, bg=NA,
           legend=c('Survey data, HIV-', 'Model, HIV-'),
           col=c('#000000', '#0066ff'), lty=1, pch=c(16,NA))
    title(xlab="Months since last birth", line=1.50)
    title(ylab="Breastfeeding, %", line=1.75)
    
    ## 2. Fit to HIV+ women, with data
    par(las=1, mar=c(2.5,3.0,1.5,1.0), cex=1)
    plot(c(0,36), c(0,1), cex=0, ann=FALSE, axes=FALSE, yaxs='i')
    # for (k in 1:nrow(draws)) {lines(months, draws[k,post.ind[,s,2]], col=cset[k])}
    polygon(c(months, rev(months)), c(fit.summary[summ.ind[2,s,],4], rev(fit.summary[summ.ind[2,s,],8])), col='#cc336630', lty=0)
    lines(months, draws[post.max,post.ind[,s,2]], col='#cc3366')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=seq(0,36,6), labels=seq(0,36,6))
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    arrows(months[data.pos$age], data.pos$lower, months[data.pos$age], data.pos$upper, code=3, angle=90, length=0.0125, col='#000000', xpd=1)
    points(months[data.pos$age], data.pos$value, pch=16, col='#000000', xpd=1)
    legend('bottomleft', inset=c(0.025,0), box.lty=0, bg=NA,
           legend=c('Survey data, HIV+', 'Model, HIV+'),
           col=c('#000000', '#cc3366'), lty=1, pch=c(16,NA))
    title(xlab="Months since last birth", line=1.50)
    title(ylab="Breastfeeding, %", line=1.75)
    title(main=sprintf("%s %s %s", abstract$CountryName[sid], abstract$SurveyYearLabel[sid], abstract$SurveyType[sid]))
    
    ## 3. Fit to HIV- and HIV+ women, don't show data
    par(las=1, mar=c(2.5,3.0,1.5,1.0), cex=1)
    plot(c(0,36), c(0,1), cex=0, ann=FALSE, axes=FALSE, yaxs='i')
    polygon(c(months, rev(months)), c(fit.summary[summ.ind[1,s,],4], rev(fit.summary[summ.ind[1,s,],8])), col='#0066ff30', lty=0)
    polygon(c(months, rev(months)), c(fit.summary[summ.ind[2,s,],4], rev(fit.summary[summ.ind[2,s,],8])), col='#cc336630', lty=0)
    lines(months, draws[post.max,post.ind[,s,1]], col='#0066ff')
    lines(months, draws[post.max,post.ind[,s,2]], col='#cc3366')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=seq(0,36,6), labels=seq(0,36,6))
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    legend('bottomleft', inset=c(0.025,0), box.lty=0, bg=NA,
            legend=c('Model, HIV-', 'Model, HIV+'),
            col=c('#0066ff', '#cc3366'), lty=1, pch=c(NA,NA))
    title(xlab="Months since last birth", line=1.50)
    title(ylab="Breastfeeding, %", line=1.75)
    
    dev.off()
  }
}

## Create a matrix of posterior mode breastfeeding model estimates. Each column
## corresponds to a survey, each row to a time since last birth. Setting
## hiv.status=1 will return estimates for HIV-negative women, hiv.status=2 for
## HIV-positive women. Results are formatted as Spectrum inputs that specify
## the probability 0 <= x <= 100 of not breastfeeding.
bf.matrix.hiv = function(stan.fit, data.flat, hiv.status) {
  ## We need to get survey ids this way so that they match the ordering in the
  ## stan fit, otherwise we will probably associate Stan's posterior model fits
  ## with the wrong surveys.
  abstract = data.flat[!duplicated(data.flat[,c('SubregionName', 'CountryName', 'SurveyYear', 'SurveyId')]),]
  
  months = cut(seq(1, 35, 2), seq(0,36,2), include.lowest=TRUE, right=FALSE)
  surveyIds = abstract$SurveyId
  n.survey = length(surveyIds)
  
  draws = as.data.frame(stan.fit)
  post.max = which.max(draws$lp__)
  post.ind = array(grep("propbf", colnames(draws)), dim=c(length(months), n.survey, 2))
  
  bf.est = matrix(nrow=length(months), ncol=nrow(abstract))
  colnames(bf.est) = surveyIds
  rownames(bf.est) = months
  
  for (s in 1:n.survey) {
    bf.est[,s] = 100 * (1.0 - unlist(draws[post.max, post.ind[,s,hiv.status]]))
  }
  
  return(bf.est)
}

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (!exists("data.flat")) {
  data.flat = load.rdhs()
  data.flat = data.flat[order(data.flat$hiv, data.flat$CountryName, data.flat$SurveyId, data.flat$age),]
  save('data.flat', file="bf-data-ws.RData") # save survey data so we can skip lengthy data download and processing later
}

post = bf.fit(data.flat)              # "post" is a misnomer since we're use MLE instead of MAP estimation
save("post", file="bf-post-ws.RData") # save model fit for later analyses

plot.model.fit(post, data.flat, "bf-post-draw.tiff")

bf.hiv = bf.matrix.hiv(post, data.flat, hiv.status=2)
write.csv(bf.hiv, "bf-post-hiv.csv")

## Write maximum likelihood parameter estimates to CSV files. These can be
## copied into Spectrum's AMModData InfantFeedingModel tab to use breastfeeding
## model estimates
mle.pars = param.tables(post, data.flat)
write.csv(mle.pars$country, "bf-pars-ceff.csv") # country effects
write.csv(mle.pars$region,  "bf-pars-reff.csv") # regional effects
