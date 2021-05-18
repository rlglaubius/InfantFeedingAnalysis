## Robert Glaubius
## Avenir Health
## 2021-05-18
##
## bf-support.R support code for bf-model.R

library(haven)
library(survey)

options(survey.lonely.psu = "remove")

set_na = function(x, na_codes = 9){x[x %in% na_codes] = NA; x}

aggr.bf.data = function(dat) {
  ## Scale weights to sum to the sample size
  obs = which(is.finite(dat$weight) & is.finite(dat$breastfeeding) & is.finite(dat$hivstatus) & is.finite(dat$monthscat))
  dat$weight = dat$weight * length(obs) / sum(dat$weight[obs])
  
  des = svydesign(ids=~v001, strata=~strata, data=subset(dat, !is.na(weight)), weights=~weight, nest=TRUE)
  
  ## Calculate the survey-weighted numbers of women who report they are/are not currently breastfeeding by HIV status
  tab = as.data.table(svytable(~breastfeeding+monthscat+hivstatus, design=des))
  val = reshape(tab, v.names='N', timevar='breastfeeding', idvar=c('monthscat', 'hivstatus'), direction='wide')
  val$monthscat = as.integer(val$monthscat)
  colnames(val) = c('age', 'hiv', 'not.bf', 'yes.bf')
  
  ## Calculate breastfeeding probabilities and confidence intervals
  est = svyby(~breastfeeding, by=~monthscat+hivstatus, design=des, svyciprop, vartype='ci')
  est$hivstatus = as.character(est$hivstatus)
  colnames(est) = c('age', 'hiv', 'value', 'lower', 'upper')
  
  num = svyby(~breastfeeding, by=~monthscat+hivstatus, design=des, unwtd.count)
  num$se = NULL
  colnames(num) = c("age", "hiv", "N")

  wgt = dplyr::left_join(val, est, by=c('age', 'hiv'))
  dat = dplyr::left_join(wgt, num, by=c('age', 'hiv'))
  return(dat)
}

load.rdhs = function() {
  surveys = as.data.table(dhs_surveys(surveyCharacteristicIds = 23))
  surveys = subset(surveys, RegionName == "Sub-Saharan Africa")
  
  ## Liberia 2019-20 DHS (LB2019DHS) did not include HIV testing, but has
  ## surveyCharacteristicIds = 23. We need to skip the survey to avoid errors.
  surveys = surveys[!SurveyId %in% c("BJ2006DHS", "BJ2012DHS", "DR2002DHS", "LB2019DHS", "ML2001DHS", "ZM2002DHS")]

  ird = as.data.table(dhs_datasets(surveyIds = surveys$SurveyId, fileType = "IR", fileFormat = "flat"))
  ard = as.data.table(dhs_datasets(surveyIds = surveys$SurveyId, fileType = "AR", fileFormat = "flat"))
  
  ird$path = unlist(get_datasets(ird$FileName))
  ard$path = unlist(get_datasets(ard$FileName))
  
  survid = surveys$SurveyId[1]
  
  ir = readRDS(ird[SurveyId == survid]$path)
  ar = readRDS(ard[SurveyId == survid]$path)
  
  ## Calculate breasfeeding by HIV status
  data.list = list()
  
  for(survid in ird$SurveyId) {
    cat(sprintf("%s...", survid))
    ir = readRDS(ird[SurveyId == survid]$path)
    
    if(survid == "CI2005AIS"){
      ## For Cote d'Ivoire 2005 AIS, individuals are uniquely identified by
      ## four variables {cluster, structure, household, line}.
      ir$v002 = 100L*ir$sstruct + ir$v002
    } else if (survid == "UG2011AIS") {
      ## Uganda 2011 AIS omits some standard DHS questions. We approximate the DHS variables of interest here.
      ## Women were asked about current breastfeeding of last born child if still alive (Q307, response in S307)
      ir$v222 = ir$v008 - ir$s212c # NOTE: Calculate months since last birth from AIS date of birth of last born child
      ir$v404 = ir$s307            # NOTE: AIS question pertains to last born child, DHS to current breastfeeding of any child
      ir$v404[is.na(ir$s307)] = 0  # We assume not currently breastfeeding if not currently breastfeeding a living last-born child
      ir$b5_01 = with(ir, ifelse(is.finite(s307), 1, NA)) # Mark lastborn child as alive if current breastfeeding asked, missing otherwise
    } else if (survid == "TZ2012AIS") {
      ## Tanzania 2012 AIS does not include v222 (interval between last birth and interview). Birth histories appear to be ordered so that b3_01 is the date of birth of the youngest child.
      ir$v222 = ir$v008 - ir$b3_01
    }
    
    if (survid %in% c("GA2012DHS", "GM2013DHS", "SN2010DHS")) {
      ir$strata = ir$v022
    } else {
      ir$strata = ir$v023
    }

    ## Select women only if combined M/R IR dataset
    if(exists("aidsex", ir)) {
      ir = subset(ir, aidsex == 2)
    }
    
    ## Check for breastfeeding variable and time since birth variable
    ## v404: currently breastfeeding
    ## v222: Last birth to interview (months)
    if(is.null(ir$v404)) {cat("skipped; v404 missing\n"); next}
    if(is.null(ir$v222)) {cat("skipped; v222 missing\n"); next}

    ar = readRDS(ard[SurveyId == survid]$path)
    if(survid == "CI2005AIS")
      ar$hivnumb = 100L*ar$hivstruct + ar$hivnumb
    names(ar) = sub("hivnumber", "hivnumb", names(ar)) # For VNAR51 dataset...
    ar[, c("hivclust", "hivnumb", "hivline")] = lapply(ar[, c("hivclust", "hivnumb", "hivline")], as.integer)
    
    ## Replace results for Zambia 2013-14 with confirmed test results
    if(survid == "ZM2013DHS")
      ar$hiv03 = ar$shiv51
    
    ## Remove labels from v001, v002 and v003 to suppress warnings during left_join
    attr(ir$v001, "label") = NULL
    attr(ir$v002, "label") = NULL
    attr(ir$v003, "label") = NULL

    dat = dplyr::left_join(ir, ar, by=c("v001" = "hivclust", "v002" = "hivnumb", "v003" = "hivline"))
    dat$weight = dat$hiv05 * 1e-6
    # dat$monthscat = as.integer(cut(ir$v222, seq(0,36,2), include.lowest=TRUE, right=FALSE)) # consistent with Spectrum
    dat$monthscat = as.integer(cut(ir$v222, seq(0,36,2), include.lowest=TRUE, right=TRUE)) # consistent with Sasi's PHIA analysis
    dat$hivstatus = factor(set_na(dat$hiv03, 4:9) > 0, c(FALSE, TRUE), c("negative", "positive"))
    dat$breastfeeding = set_na(dat$v404, 9)
    
    if(all(is.na(dat$breastfeeding))) {cat("skipped; v404 contains no data\n"); next}
    
    ## b5_01==1 restricts to women whose last-born child is still alive
    data.list[[survid]] = aggr.bf.data(subset(dat, b5_01==1))
    data.list[[survid]]$age = as.integer(data.list[[survid]]$age)
    
    data.list[[survid]]$SurveyId        = survid
    data.list[[survid]]$CountryName     = surveys[SurveyId == survid]$CountryName
    data.list[[survid]]$SurveyYear      = surveys[SurveyId == survid]$SurveyYear
    data.list[[survid]]$SurveyYearLabel = surveys[SurveyId == survid]$SurveyYearLabel
    data.list[[survid]]$SubregionName   = surveys[SurveyId == survid]$SubregionName
    data.list[[survid]]$SurveyType      = surveys[SurveyId == survid]$SurveyType

    cat("done\n")
  }
  
  data.flat = do.call(rbind, data.list)
  return(data.flat[order(data.flat$hiv, data.flat$CountryName, data.flat$SurveyId, data.flat$age),])  
}

param.tables = function(stan.fit, data.flat) {
  abstract = data.flat[!duplicated(data.flat[,c('SubregionName', 'CountryName', 'SurveyYear', 'SurveyId')]),]
  
  username = Sys.info()['user']
  code.map = read_excel(sprintf('C:/Users/%s/Dropbox (Avenir Health)/My Files/R/dhs-iso3166-map.xlsx', username), sheet='CodeMap')
  
  draws = as.data.frame(stan.fit)
  mle.ind = which.max(draws$lp__)
  
  country.names = levels(factor(abstract$CountryName))
  country.abbrs = substr(abstract$SurveyId[match(country.names, abstract$CountryName)], 1, 2)
  country.codes = sapply(country.abbrs, function(code) {code.map$ISO.Numeric[which(code.map$DHS.Code==code)]})
  country.region = abstract$SubregionName[match(country.names, abstract$CountryName)]
  
  region.names = levels(factor(abstract$SubregionName))
  
  ## data frame of country-level parameter values
  cpar.df = data.frame(
    country = country.names,
    isocode = country.codes,
    region  = country.region,
    bgn = unlist(draws[mle.ind, grep('ceff_bgn', names(draws))]), # initial BF probability
    med = unlist(draws[mle.ind, grep('ceff_med', names(draws))]), # BF duration median among mothers who initially BF
    shp = unlist(draws[mle.ind, grep('ceff_shp', names(draws))])) # BF duration shape among mothers who initially BF
  
  rpar.df = data.frame(
    region = region.names,
    bgn.time     = unlist(draws[mle.ind, grep('base_slope_bgn', names(draws))]), # initial BF probability: time effect
    bgn.hiv      = unlist(draws[mle.ind, grep('base_bgn_hiv',   names(draws))]), # initial BF probability: HIV effect
    bgn.hiv.time = unlist(draws[mle.ind, grep('base_bgn_arv',   names(draws))]), # initial BF probability: interaction between HIV status and time
    med.time     = unlist(draws[mle.ind, grep('base_slope_med', names(draws))]), # median BF: time effect
    med.hiv      = unlist(draws[mle.ind, grep('base_med_hiv',   names(draws))]), # median BF: HIV effect
    med.hiv.time = unlist(draws[mle.ind, grep('base_med_arv',   names(draws))]), # median BF: interaction between HIV status and time
    shp.time     = unlist(draws[mle.ind, grep('base_slope_shp', names(draws))])) # BF shape: time effect
  
  return(list(
    country = cpar.df,
    region  = rpar.df))
}

