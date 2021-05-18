## Robert Glaubius
## Avenir Health
## 2021-05-18
##
## gen-figure.R script to generate figure showing estimated breastfeeding
## patterns over time, by region and maternal HIV status. This expects that
## survey data have already been downloaded and the breastfeeding model has been
## fit, with corresponding data available in bf-data-ws.RData and
## bf-post-ws.RData, respectively. If those are not present, you must run
## source("bf-model.R") to generate them before running this script.

library(boot)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(rstan)

## Path to stored model fit RData files 
path.rdata = "."
load(sprintf("%s/bf-data-ws.RData", path.rdata)) # load saved survey data
load(sprintf("%s/bf-post-ws.RData", path.rdata)) # load saved model fit

annotate.posterior = function(model.data, model.post) {
  longify.ceff = function(par.cols, country.names, country.alpha, country.region) {
    pars.wide = cbind.data.frame(country = country.names, alpha = country.alpha, region = country.region, t(par.cols))
    pars.long = reshape(pars.wide,
                        varying = rownames(par.cols),
                        v.names = "value",
                        timevar = "sample",
                        times   = 1:nrow(par.cols),
                        direction = "long")
    rownames(pars.long) = NULL
    pars.long$id = NULL
    return(pars.long)
  }

  longify.reff = function(par.cols, region.names) {
    pars.wide = cbind.data.frame(region = region.names, t(par.cols))
    pars.long = reshape(pars.wide,
                        varying = rownames(par.cols),
                        v.names = "value",
                        timevar = "sample",
                        times   = 1:nrow(par.cols),
                        direction = "long")
    rownames(pars.long) = NULL
    pars.long$id = NULL
    return(pars.long)
  }

  abstract = model.data[!duplicated(model.data[,c('SubregionName', 'CountryName', 'SurveyYear', 'SurveyId')]),]
  code.map = read_excel('dhs-iso3166-map.xlsx', sheet='CodeMap')

  country.names = levels(factor(abstract$CountryName))
  country.abbrs = substr(abstract$SurveyId[match(country.names, abstract$CountryName)], 1, 2)
  country.alpha = sapply(country.abbrs, function(code) {code.map$ISO.Alpha.3[which(code.map$DHS.Code==code)]})
  country.region = abstract$SubregionName[match(country.names, abstract$CountryName)]

  region.names = levels(factor(abstract$SubregionName))

  ## Label country effects
  ceff = list(
    bgn = longify.ceff(model.post[,grep("ceff_bgn", names(model.post))], country.names, country.alpha, country.region),
    med = longify.ceff(model.post[,grep("ceff_med", names(model.post))], country.names, country.alpha, country.region),
    shp = longify.ceff(model.post[,grep("ceff_shp", names(model.post))], country.names, country.alpha, country.region))

  ## Label regional effects
  reff = list(
    bgn_time     = longify.reff(model.post[,grep("base_slope_bgn", names(model.post))], region.names),
    bgn_hiv      = longify.reff(model.post[,grep("base_bgn_hiv",   names(model.post))], region.names),
    bgn_hiv_time = longify.reff(model.post[,grep("base_bgn_arv",   names(model.post))], region.names),
    med_time     = longify.reff(model.post[,grep("base_slope_med", names(model.post))], region.names),
    med_hiv      = longify.reff(model.post[,grep("base_med_hiv",   names(model.post))], region.names),
    med_hiv_time = longify.reff(model.post[,grep("base_med_arv",   names(model.post))], region.names),
    shp_time     = longify.reff(model.post[,grep("base_slope_shp", names(model.post))], region.names))

  return(list(ceff = dplyr::bind_rows(ceff, .id="param"),
              reff = dplyr::bind_rows(reff, .id="param")))
}

post.flat = as.data.frame(post)
post.flat.long = annotate.posterior(data.flat, post.flat)

mle.first = which.max(post.flat$lp__)

effects.regional = post.flat.long$reff
effects.regional$region = factor(plyr::mapvalues(effects.regional$region,
                                                 from = c("Eastern Africa", "Middle Africa",  "Southern Africa", "Western Africa"),
                                                 to   = c("Eastern Africa", "Central Africa", "Southern Africa", "Western Africa")))

effects.national = post.flat.long$ceff
effects.national$region = factor(plyr::mapvalues(effects.national$region,
                                                 from = c("Eastern Africa", "Middle Africa",  "Southern Africa", "Western Africa"),
                                                 to   = c("Eastern Africa", "Central Africa", "Southern Africa", "Western Africa")))

average.regional = plyr::ddply(effects.national, .(region, param, sample), function(df) {
  data.frame(region  = df$region[1],
             param   = df$param[1],
             sample  = df$sample[1],
             value   = mean(df$value))
})

regional.params = rbind(average.regional, effects.regional)
regional.params.wide = reshape(regional.params,
                               timevar = "param",
                               idvar = c("region", "sample"),
                               direction = "wide")
colnames(regional.params.wide) = gsub("value.", "", colnames(regional.params.wide))

regional.pattern = plyr::ddply(regional.params.wide, .(sample, region), function(df) {
  dplyr::bind_rows(lapply(c(2015, 2010, 2005), function(year) {
    t.diff = year - 2010
    bgn.neg = with(df, inv.logit(logit(bgn) + bgn_time * t.diff))
    med.neg = with(df, exp(med + med_time * t.diff))

    bgn.pos = with(df, inv.logit(logit(bgn) + (bgn_time + bgn_hiv_time) * t.diff + bgn_hiv))
    med.pos = with(df, exp(med + (med_time + med_hiv_time) * t.diff + med_hiv))
    shp = with(df, exp(shp + shp_time * t.diff))

    x = 0:36
    y.neg = bgn.neg * (1.0 - 1.0 / (1.0 + (x / med.neg)^(-shp)))
    y.pos = bgn.pos * (1.0 - 1.0 / (1.0 + (x / med.pos)^(-shp)))

    rbind(data.frame(sample = df$sample[1],
                     region = df$region[1],
                     year  = year,
                     hiv   = "HIV-",
                     age   = x,
                     value = y.neg),
          data.frame(sample = df$sample[1],
                     region = df$region[1],
                     year  = year,
                     hiv   = "HIV+",
                     age   = x,
                     value = y.pos))
  }))
})

regional.stats = plyr::ddply(regional.pattern, .(region, year, hiv, age), function(df) {
  value = df$value[df$sample==mle.first]
  cred = quantile(df$value, c(0.025, 0.975))

  data.frame(
    region  = df$region[1],
    year    = df$year[1],
    age     = df$age[1],
    hiv     = df$hiv[1],
    value   = value,
    lower   = cred[1],
    upper   = cred[2])
})
regional.stats$region = factor(regional.stats$region, levels=c("Western Africa", "Central Africa", "Eastern Africa", "Southern Africa"))
regional.stats$year = factor(regional.stats$year)

cset = brewer.pal(8, "Set1")[1:3]
ggplot(regional.stats, aes(x=age, y=100*value, ymin=100*lower, ymax=100*upper, group=year)) +
  geom_ribbon(aes(fill=year), alpha=0.4) +
  geom_line(aes(color=year)) +
  facet_grid(hiv~region) +
  ylim(0,100) +
  ylab("Mothers currently breastfeeding, %") +
  scale_x_continuous(name="Child's age in months", breaks=seq(0,36,6)) +
  scale_fill_manual(values=cset) +
  scale_color_manual(values=cset) +
  theme_bw() +
  theme(legend.position="right",
        legend.margin=margin(l=0, r=0, unit="cm"),
        plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
        panel.border = element_rect(fill=NA, color="#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size=0.4),
        strip.background = element_rect(fill="grey85", color="#000000"),
        text = element_text(size=10),
        axis.title = element_text(size=rel(1.0)),
        axis.text = element_text(color="#000000", size=rel(0.8)))
ggsave("figure-breastfeeding.tiff", compression="lzw", dpi=600, units="mm", width=180, height=60)

