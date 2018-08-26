######################################################################################################################
# Author: Alex Keil
# Program: workshop_code.sas
# Date: 20180826
# Project: ISEE 2018 workshop: Causal inference foundations and applications in environmental health sciences
# Tasks:
# Data in: miners, coalplant
# Description: Carry out parametric g-formula analysis and strucgtural nested model analysis of data simulated to
#   emulate an occupational study of uranium miners looking at the effects of occupational exposures to radon
#   and all cause mortality, as well as a population study of air pollutants that arise from coal-fired power
#   plants
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

# important packages
# install with: 
# install.packages(c('ggplot2', 'readr', 'dplyr', 'survival', 'tibble', 'boot'), repos = 'https://cran.mtu.edu')
library(ggplot2)
library(readr)
library(dplyr)
library(survival)
library(tibble)
library(boot)
# reading in data
coalplant = read_csv("~/2018_ISEE_causal/analyses/data/coalplant.csv")
miners = read_csv("~/2018_ISEE_causal/analyses/data/miners.csv")
outpath = "~/2018_ISEE_causal/analyses/output/"

################################################################################
#  Coal plant example - parametric g-formula with point exposure
################################################################################
# MDI (mental development index) appears to be roughly normally distributed
#  so we will model using a linear model

mdimod = glm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=coalplant)

# Say we know based on local emissions data and discussions with atmospheric
#  chemists that 91% of As, 96% of Be, and 45% of Cd ambient levels come
#  from a local coal-fired power plant


gformula_means <- function(ymodel, data, ...){
   # function to calculate mean MDI under a given intervention
   require('dplyr')
   # if only 'ymodel' and 'data' are supplied, then natural course is fit
  postinterventionX <- mutate(data,...)
  mean(predict(ymodel, newdata=postinterventionX))
}

#point estimates for mean MDI under each intervention
nc = gformula_means(mdimod, data=coalplant)
no_coal = gformula_means(mdimod, data=coalplant, 
                         cd=cd*(1-0.5), as=as*(1-0.96), be=be*(1-0.91))
no_cd = gformula_means(mdimod, data=coalplant, cd=0)
no_as = gformula_means(mdimod, data=coalplant, as=0)
no_be = gformula_means(mdimod, data=coalplant, be=0)
# expected population mean MDI
print(c(nc=nc, no_coal=no_coal, no_cd=no_cd, no_as=no_as, no_be=no_be), 4)

# bootstrap to get confidence intervals on effects (boot function also gives point estimates from original data)
gformula_meandiffs <- function(data, index){
  mcsample = data[index,]
  bootmod = glm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=mcsample)
  nc = gformula_means(bootmod, data=mcsample)
  no_coal = gformula_means(bootmod, data=mcsample, cd=cd*(1-0.5), as=as*(1-0.96), be=be*(1-0.91))
  no_cd = gformula_means(bootmod, data=mcsample, cd=0)
  no_as = gformula_means(bootmod, data=mcsample, as=0)
  no_be = gformula_means(bootmod, data=mcsample, be=0)
  c(nc_mdi=nc, shutdown=no_coal-nc, attrmdi_cd=no_cd-nc, attrmdi_as=no_as-nc, attrmdi_be=no_be-nc)
}

set.seed(12321)
bootsamples = boot(data=coalplant, statistic=gformula_meandiffs, R=500)

# effect estimates with confidence intervals
se = apply(bootsamples$t, 2, sd)
print(cbind(estimate=bootsamples$t0, lowerCI=bootsamples$t0-1.96*se, upperCI=bootsamples$t0+1.96*se), 3)

################################################################################
#  Uranium miners examples
################################################################################



################################################################################
#
#
# g-computation/g-formula
#
#
################################################################################

get_coefs <- function(data){
  # pooled linear model for log(exposure), while at work
  data$lograd = with(data, ifelse(atwork==1, log(rad), NA))
  mod_e = glm(lograd ~ outtime + cum_rad_lag1 + smoker, data=filter(data, atwork==1))
  # pooled logistic model for leaving work (at the end of the year), while at work
  mod_w = glm(leavework ~ outtime + I(outtime^2) + rad_lag1 + cum_rad_lag1 + smoker, data=filter(data, atwork==1 | leavework==1), family=binomial(link=logit))
  # pooled logistic model for death
  mod_d = glm(dead ~ outtime + atwork + cum_rad + I(cum_rad*cum_rad) + smoker + smoker*cum_rad, data=data, family=binomial(link=logit))
  # when troubleshooting code, it is helpful to have a time-only model (which should fit the natural course well)
  #mod_e = glm(lograd ~ outtime + I(outtime^2), data=filter(data, atwork==1))
  #mod_w = glm(leavework ~ outtime + I(outtime^2), data=filter(data, atwork==1 | leavework==1), family=binomial(link=logit))
  #mod_d = glm(dead ~ outtime + I(outtime^2) + I(outtime^3), data=data, family=binomial(link=logit))
  list(mod_e=mod_e, mod_w=mod_w, mod_d=mod_d)
}

gformula_risks <- function(intervention=NULL, pseudo_cohort_bl, 
                           endtime, mod_e, mod_w, mod_d, seed=NULL){
  # esimate intervention specific risks using a Monte Carlo algorithm
  # this function assumes ordering by ID, TIME
  require(dplyr)
  N <- dim(pseudo_cohort_bl)[1]
  ncols <- dim(pseudo_cohort_bl)[2]
  # create container to hold simulated results
  pseudo_cohort <- slice(pseudo_cohort_bl,rep(1:N, each=endtime))
  pseudo_cohort$gfid <- rep(seq(1, N, 1), each=endtime) # new id for each pseudo individual
  
  # initialize values
  pseudo_cohort$intime <- rep(seq(0, (endtime-1), 1), times=N)
  pseudo_cohort$outtime <- pseudo_cohort$intime+1
  pseudo_cohort$atwork <- 1
  pseudo_cohort$durwork <- 1
  pseudo_cohort$leavework <- 0
  pseudo_cohort$rad <- 0
  pseudo_cohort$rad_lag1 <- 0
  pseudo_cohort$cum_rad_lag1 <- 0
  pseudo_cohort$durwork_lag1 <- 0
  pseudo_cohort$dead <- 0
  pseudo_cohort$cum_hazard <- 0 # new cumulative hazard
  pseudo_cohort$cum_incidence <- 0 # new cumulative incidence
  # error term variance for exposure model
  #var_e = sum(mod_e$residuals^2)/mod_e$df.residual # default R approach (unbiased SD - better approach)
  var_e = mean(mod_e$residuals^2) # to match with SAS version (unbiased variance)
  set.seed(seed)
  # limits on log exposure (keep simulated values in range of observed)
  max_e = 7.321276
  
  # simulate time-varying values
  # R is slow with looping, so this simulation is vectorized as much as possible
  for(t in seq(1, endtime, 1)){
    # index that keeps track of time
    idx = which(pseudo_cohort$outtime == t)
    #update all values of lagged variables (easy to forget!)
    if(t > 1){
      pseudo_cohort[idx,'cum_rad_lag1'] = pseudo_cohort[idx-1,'cum_rad']
      pseudo_cohort[idx,'rad_lag1'] = pseudo_cohort[idx-1,'rad']
    }    
    ############ 
    # leaving work (time varying confounder)
    ############
    if(t==1) widx = 1:length(idx) # everyone is at work at first time point
    if(t > 1){
       widx = which(pseudo_cohort[idx-1,'atwork']==1)
       # index keeping track of which workers still work
       pseudo_cohort[idx[widx],'leavework'] <- rbinom(n=length(widx), size=1, 
                                                      prob=predict(mod_w, newdata = pseudo_cohort[idx[widx],], 
                                                                   type = 'response'))
       # if worker didn't leave, then assume stay at work
       pseudo_cohort[idx,'atwork'] <- (pseudo_cohort[idx,'leavework']==0 & pseudo_cohort[idx-1,'atwork']==1)
       # update at work index to account for workers who left
       widx = which(pseudo_cohort[idx,'atwork']==1)
       pseudo_cohort[idx,'durwork'] <- pseudo_cohort[idx-1,'durwork'] + pseudo_cohort[idx,'atwork']
       pseudo_cohort[idx,'durwork_lag1'] <- pseudo_cohort[idx-1,'durwork']
    }
    ############
    # exposure/interventions
    ############
    # exposure is assumed log-normally distributed (=predicted value plus draw from the residuals)
    meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[widx],])
    logr = meanlogr + rnorm(n=length(widx), mean=0, sd=sqrt(var_e))
    pseudo_cohort[idx[widx],'rad'] <- exp(pmin(log(max_e), logr))
    # exposure under interventions
    if(typeof(intervention) != "NULL"){
      if(is.numeric(intervention)){
        # static, deterministic intervention, e.g. set exposure = 0
        pseudo_cohort[idx,'rad'] <- intervention
      } else{
      if(is.character(intervention)){
        # dynamic interventions are defined by character strings: e.g. 'rad>0.3' means we intervene to prevent radon exposure from being > 0.3
        #  this line creates an index of rows at the current time which are in violation of the intervention
        #  exposure will be redrawn for these observations until they are in compliance.
        #  This is equivalent to drawing from a truncated log-normal distribution
        viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
        while(length(viol_idx)>0){
          # accept/rejection sampling to get draws from truncated log-normal
          # this is how we (I) assume an intervention would be implemented, but
          # we could imagine other ways (e.g. a hard cap on exposure)
          meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[viol_idx],])
          lograd = meanlogr + rnorm(n=length(viol_idx), mean=0, sd=sqrt(var_e))
          pseudo_cohort[idx[viol_idx],'rad'] <- exp(lograd)
          # check whether any are still in violation of the distribution, if so, repeat loop
          viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
        }}} # end dynamic intervention
    } # end all interventions on exposure
    if(t > 1){
      pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx-1,'cum_rad'] + pseudo_cohort[idx,'rad']
    } else pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx,'rad']
    ############
    # death (simulate discrete hazard as in Taubman et al 2009, rather than 1/0 as in Keil et al 2014)
    #  see note at bottom of program about late entry - in the case of late entry
    #  the typical approach is to simulate actual death (1/0 variable) and then
    #  estimate survival in the pseudo-population using a survival curve estimator
    #  that can account for late entry
    ############
    pseudo_cohort[idx,'dead'] = predict(mod_d, newdata = pseudo_cohort[idx,], type = 'response')
    if(t > 1){
      # product limit estimator of cumulative incidence (note this is applied on the individual basis)
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx-1,'cum_incidence'] + 
        (1-pseudo_cohort[idx-1,'cum_incidence'])*pseudo_cohort[idx,'dead']
    } else{
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx,'dead']
    }
    #alternative estimator of the cumulative incidence: Breslow estimator
    #if(t > 1){
    #  # kaplain-meier estimator of cumulative incidence (note this is applied on the individual basis)
    #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx-1,'cum_hazard'] + pseudo_cohort[idx,'dead']
    #} else{
    #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx,'dead']
    #}
    #  pseudo_cohort[idx,'cum_incidence'] = 1-exp(-pseudo_cohort[idx,'cum_hazard'])
    # could also use Breslow estimator, but need aalen-johansen estimator if there are competing risks
    # the 'cum_incidence' variable is an estimate of the individual risk. We then
    # average that over the population below to get the marginal risk
  } # end time loop
  pseudo_cohort
}# end function

# observed data
baseline_data <- filter(miners, intime==0)
N = dim(baseline_data)[1]


# large sample from observed baseline data
mc_iter = 20000
set.seed(12325)
mcdata <- slice(baseline_data, sample(1:N, size = mc_iter, replace=TRUE))

# example of getting the cumulative incidence curves
mods = get_coefs(data=miners)
nc = gformula_risks(intervention=NULL, pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
unex = gformula_risks(intervention=0, pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
exp3 = gformula_risks(intervention='rad>0.3', pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
exp1 = gformula_risks(intervention='rad>0.1', pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)

# comparing natural course and observed, cumulative incidence
# we get population risk/cumulative incidence by the mean of the individual risks
obscifit <- survfit(Surv(intime, outtime, dead)~1, data=miners)
gfcinc <- with(nc, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciunex <- with(unex, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciexp3 <- with(exp3, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciexp1 <- with(exp1, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
obsnc <- tibble(time=obscifit$time, ci=1-obscifit$surv)

# time ratio for comparison with SNM
# total expected years of life in natural course versus intervention to limit exposure to be below 0.1
# with no censoring this calculation is easy! The total expected number of years lived is just
# the survival function summed over all the person-time.
elifenc = with(nc, sum(1-cum_incidence))
elifeexp1 = with(exp1, sum(1-cum_incidence))
elifeexp1/elifenc # [1] 1.476379
# recall g-estimate was 1.36 (1.09, 1.67)


ggplot(aes(x=time, y=ci), data=obsnc) + 
  geom_line(aes(color="Observed")) +
  geom_line(aes(x=time, y=ci, color="Natural course"), data=gfcinc) +
  geom_line(aes(x=time, y=ci, color="Unexposed"), data=gfciunex) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.3"), data=gfciexp3) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.1"), data=gfciexp1) +
  scale_y_continuous(name="Cumulative incidence") +
  scale_x_continuous(name="Time", limits=c(0,20)) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(legend.position = c(0.99,0.01),legend.justification = c(1,0))+
  ggtitle(expression(paste("G-formula estimates of cumulative incidence")))
ggsave(paste0(outpath, 'r_gformulaCI.png'), width=5, height=4)
  


# now we turn the g-formula algorithm to get risk difference into a function in order
# to get bootstrap variance
gformula_effectestimates <- function(data=baseline_data, index=NULL, fdata=miners, endtime=10, seed=NULL){
  # for bootstrap samples, when needed: first take a random sample of the 
  #  study IDs, with replacement
  # 0) preprocessing: take bootstrap sample if bootstrapping, otherwise use the observed data
    if(is.null(index)){
      # if index is NULL, then just use original data and do nothing here
    } else{
      # bootstrap sample of the baseline data using the vector 'index'
      data = slice(data, index)
      # resample from the full data according to the id variable in the baseline data
      idindex = tibble(reps = table(data$id), id = as.numeric(names(table(data$id))))
      fdata = merge(idindex, fdata, by='id')
      fdata = slice(fdata, rep(1:dim(fdata)[1], times=fdata$reps))
    }
  # 1) fit models to full data/bootstrap sample of the full data
    mods = get_coefs(data=fdata)
  # 2) take large MC sample from baseline data
    mc_iter = 20000
    set.seed(seed)
    mcdata <- slice(data, sample(1:N, size = mc_iter, replace=TRUE))
  # 3) simulate probability distributions in large MC sample using model fits from step 1
    nc = gformula_risks(intervention=NULL,        pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    unex = gformula_risks(intervention=0,         pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    exp3 = gformula_risks(intervention='rad>0.3', pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    exp1 = gformula_risks(intervention='rad>0.1', pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    # note: I change the seed for each of these to emulate SAS approach - this is not required
  # 4) summarize data
    ci_nc = as.numeric(with(nc, tapply(cum_incidence, outtime, mean)))
    ci_unex = as.numeric(with(unex, tapply(cum_incidence, outtime, mean)))
    ci_exp3 = as.numeric(with(exp3, tapply(cum_incidence, outtime, mean)))
    ci_exp1 = as.numeric(with(exp1, tapply(cum_incidence, outtime, mean)))
   # return a vector of risk differences and log-risk ratios at time t=endtime
   # could also calculate for all times
     c(
       rd_exp3_10=c(ci_exp3-ci_nc)[endtime],
       rd_exp1_10=c(ci_exp1-ci_nc)[endtime],
       rd_unex_10=c(ci_unex-ci_nc)[endtime],
       lrr_exp3_10=log(c(ci_exp3/ci_nc)[endtime]),
       lrr_exp1_10=log(c(ci_exp1/ci_nc)[endtime]),
       lrr_unex_10=log(c(ci_unex/ci_nc)[endtime])
     )
}

# point estimate for risk differences at t=10 (also given by the 'boot' function below)
system.time(gf_est <- gformula_effectestimates(data=baseline_data, fdata=miners, endtime=10, seed=NULL))
# took about 3 seconds on my machine
# risk difference, log RR estimates (will vary slightly by seed value due to simulation approach)
gf_est

# getting confidence intervals: non-parametric bootstrap
# note that to get valid SE, we should use 200+ bootstrap samples
# for percentile based confidence intervals we should use 1000+
# for an approach that takes 3 seconds to run, 200 bootstrap samples = 10 minutes
#nbootsamples = 200
#system.time(boot_samples <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = nbootsamples, endtime=10))
#   user  system elapsed 
#674.059 227.130 901.453 


nbootsamples = 10
system.time(boot_samples <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = nbootsamples, endtime=10))
boot_samples$t0

# point estimate and confidence intervals
est = boot_samples$t0
lci =  boot_samples$t0 - 1.96*apply(boot_samples$t, 2, sd)
uci =  boot_samples$t0 + 1.96*apply(boot_samples$t, 2, sd)
cat(paste0(names(gf_est), ":", round(est, 3), " (", round(lci, 3), ", ", round(uci, 3), ")\n"))



# Further considerations:
# Administrative censoring: generally, we would stop follow-up of the pseudo-cohort
#   at the time of administrative censoring, but if we assume we have correct models
#   then follow-up could be, in principle, extended.
# Censoring by loss-to-follow-up: under assumption of correct model specification
#   the approach to censoring in g-computation is to simulate individuals covariate
#   and event history after the time of censorin
# Competing risks: The discrete, cause specific hazard from g-computation will be
#   consistent under correct model specification, even with competing risks. 
#   We can also "account" for competing risks by using an estimate of the cumulative
#   incidence that allows for competing events, like the Aalen-Johansen estimator.
#   If we use the Kaplan-Meier/product limit estimator with competing risks, we estimate the 
#   "conditional" risk rather than the unconditional risk. The unconditional risk
#   is typically what we want.
# Late entry: notice that the pseudo-cohort has all individuals enter at the same
#   time. This time scale is "time since hire" or "time on study" if this study
#   included individuals hired prior to the start of the study. We could late
#   enter people for a different time scale (e.g. age), but would have to calculate
#   risk differently - we would have to simulate values for death and then use an
#   estimator of the cumulative incidence that allows for late entry, such as
#   the Kaplan meier estimator
# Last: Functional form of all of the models matter. We could fit, for example,
#   machine learning or data adaptive approaches


################################################################################
# structural nested accelerated failure time model (fit with g-estimation)
# an SNM is a model for the effect of a 'blip' of treatment, conditional on being untreated in the future
# it compares the outcomes we would observe under no treatment to the outcomes we would observe
#   if treated now and then kept untreated in the future. This is an admittedly odd
#   way to think about causal effects, but it corresponds to a couple things we might be familiar with:
#   1) in time fixed settings, the only 'blip' is the time-fixed exposure so a standard (say) linear
#      model for a continuous outcome, adjusted for covariates is equivalent to a conditional 
#      structural nested model
#   2) in time-varying exposure settings, if we add all of the blip effects for "always treated"
#      we estimate the same "always treated" counterfactual as we do in other causal methods.
#      In other words, the "blip effect" is just a way of decomposing the "always" treated effect
#      
################################################################################


blipdown <- function(psi, duration, exposure){
  # blip function: the effect of a blip of exposure, holding future exposure
  #  to equal zero. Calculated for each row of the dataset.
  # log-linear SNM
  # observed "duration" of time for observation (e.g. 1 person year)
  # is multiplied by exp(X*psi), which gives the expected 
  #  contribution to the total survival time under no exposure
  pyr0 = exp(psi*exposure)*duration
  pyr0
}

t_calc <- function(id, starttime, pyr0, method=1){
  # calculate total expected lifetime under no exposure from time 'start'
  # until the end of the observed time. Calculated for each observation
  # in the dataset
  # the method described in the theory
  if(method==1) t_blip = starttime + unlist(tapply(pyr0, id, function(x) rev(cumsum(rev(x)))))
  # a slightly less efficient method that is easier in most software
  # this calculates the potential outcome under "never exposed" for all observations
  # See Hernan et al (2005) Pharmacoepidem Dr S
  if(method==2) t_blip =  rep(tapply(pyr0, id, sum), times= tapply(pyr0, id, length))
  t_blip
}

g_function_psi <- function(exp.model.formula, t_blip, data, ...){
  # calculate the 'g-function', or the test-statistic comparing observed
  #  exposure with the blipped-down time (this is a function of a particular
  #  value of psi and does not estimate psi!)
  f = as.formula(sub("~", paste("~", t_blip, "+"), x = exp.model.formula))
  mod = glm(f, data=data, ...)
  chisq <- summary(mod)$coefficients[2,3]^2
  chisq
}


################################
# g-estimation of a dose-response
################################

g_estimation <- function(psi.seq = seq(0, 1, 0.01)){
  # the main function
  chisq = numeric(length(psi.seq))
  for(i in 1:length(psi.seq)){
    miners$pyr0 <- blipdown(psi.seq[i], miners$outtime-miners$intime, miners$cum_rad)
    miners$t_blip <- t_calc(miners$id, start=miners$intime, miners$pyr0, method=1)
    #miners$t_blip <- t_calc(miners$id, start=miners$intime, miners$pyr0, method=2)
    # note that the exposure model is fit only to the person time for which individuals
    #  are at work! This leverages the "selective ignorability" condition, where
    #  we need only test exchangeability in a subset of the data in which it is 
    #  expected to hold. In this case, we necessarily fit the model to the employed
    #  person time because exposure does not occur outside of work (i.e. non-positivity)
    # Model: pooled linear model for log(exposure) as a function of time (quadratic), 
    #  cumulative prior exposure, employment duration, and smoking. In practice, we would usually fit
    #  more complex, flexible functions (e.g. splines on continuous terms, interaction terms)
    # that would reduce model form assumptions (at the cost of some extra variance).
    chisq[i] = g_function_psi(exp.model.formula='log(rad) ~ intime + I(intime^2) + cum_rad_lag1 + smoker', t_blip='t_blip', 
                              data=filter(miners, atwork==1), family=gaussian())
  }
  tibble(psi=psi.seq, chisq=chisq)
}

# The g-function is a series of 'objective function' (chi-squared values) 
#  across a grid of possible values of psi - this is a 
#  quick and dirty way to get g-estimates. We could use
#  an optimization routine (e.g. nelder-mead simplex algorithm)
#  to find the minimum, but a grid search over a set of plausible values is a 
#  more sure way to find the minimum since the objective function may have 
#  local minima - this is challenging when there are multiple psi parameters,
#  so optimization may be a better option in that case
# Here: we get chi-squared values for a grid of possible values of psi ranging
#  from -1 to -0.4 (the data are simulated, so I know the true estimate lies in this
#  range. In practice, we might want to have a larger range.
gfunction = g_estimation(psi.seq = seq(0, 1, 0.01))
min(gfunction$chisq) # should be close to zero

# Plotting the g-function. the g-estimate of psi is given by the value of psi for
# which the following hold:
#  1) the calculated potential outcomes (given by the structural model) are 
#    statistically independent from the observed exposure (given the exposure model)
#  2) the p-value for the association between the calculated potential outcomes and the
#     observed exposure is 1.0, in a test of association that is adjusted for
#     time-varying and baseline confounders
#  or 3) the coefficient for a variable with the calculated potential outcomes is zero
#     in a model for the exposure, controlling for time-varying and baseline confounders
#  These conditions are stated in order of specificity: we could establish independence
#   using (for example) log-rank tests. In this example, we use a regression model for the
#   log(exposure) [NOT cumulative exposure], adjusted for time-fixed and time-varying confounders.
#   Note that this is simply an easy (relatively) to understand approach and is not 
#   the most computationally efficient!
#   In the exposure model, we use a Wald statisic (coefficient^2/variance) ~ Chi^2(df=1) to test
#   whether the coefficient for the (guessed) potential outcomes are independent of observed
#   exposure. By conditional exchangeability, the calculated potential outcomes and the exposure
#   will be independent at the true value of psi. A confidence interval set is created by
#   inverting the Chi-squared statistic to find the points closest to the g-estimate of psi
#   at which the Chi-squared statistic is equal to the critical value of 3.84 (we could be
#   more exact and use linear interpolation between the grid-points, but an approximation will 
#   suffice here).

(g_estimate <- gfunction$psi[which.min(gfunction$chisq)]) # 0.59, true simulated value = 0.7


ggplot(data=gfunction, aes(x=psi, y=chisq)) + 
  scale_y_continuous(name=expression(chi['df=1']^2)) + 
  scale_x_continuous(name=expression(psi)) + 
  geom_step() + 
  geom_vline(xintercept=g_estimate) + 
  geom_hline(yintercept=qchisq(.95, 1)) + 
  geom_text(aes(y=qchisq(.95, 1), x=0.25, label="95% CI"), nudge_y=1) + 
  geom_text(aes(y=1, x=g_estimate, label='hat(psi)'), nudge_x=0.02, parse=TRUE) + 
  theme_classic()+
ggtitle(expression(paste("G-estimate of ", psi, " and 95% confidence intervals")))
ggsave(paste0(outpath, 'r_gfunction.png'), width=5, height=4)

# determine the confidence intervals by finding the values of psi corresponding to the 
#  closest chi-squared values to the 95th percentile of a chi-squared distribution
#  with the number of df equal to the number of coefficients tested. This is
#  'inverting the test statistic' to find the confidence set
lci = gfunction$psi[gfunction$psi<g_estimate][which.min((gfunction$chisq[gfunction$psi<g_estimate]-qchisq(.95, 1))^2)]
uci = gfunction$psi[gfunction$psi>g_estimate][which.min((gfunction$chisq[gfunction$psi>g_estimate]-qchisq(.95, 1))^2)]

# the 'time ratio' and associated confidence intervals
cat(paste0(signif(exp(g_estimate), 2), " (", signif(exp(lci), 2), ', ', signif(exp(uci), 2), ")"))
# 1.8 (1.4, 2.2)


################################
# g-estimation of policy effects
################################
g_estimation_policy <- function(psi.seq = seq(0, 1, 0.01), occlimit){
  # the main function
  chisq = numeric(length(psi.seq))
  for(i in 1:length(psi.seq)){
    miners$radlt0_1 = (miners$rad>occlimit)
    miners$pyr0p <- blipdown(psi.seq[i], miners$outtime-miners$intime, miners$radlt0_1)
    miners$t_blip <- t_calc(miners$id, start=miners$intime, miners$pyr0p, method=1)
    #miners$t_blip <- t_calc(miners$id, start=miners$intime, miners$pyr0, method=2)
    # note that the exposure model is fit only to the person time for which individuals
    #  are at work! This leverages the "selective ignorability" condition, where
    #  we need only test exchangeability in a subset of the data in which it is 
    #  expected to hold. In this case, we necessarily fit the model to the employed
    #  person time because exposure does not occur outside of work (i.e. non-positivity)
    # Model: pooled logistic model for exposure above or below the proposed limit
    #  as a function of time (quadratic), 
    #  cumulative prior exposure, employment duration, and smoking. In practice, we would usually fit
    #  more complex, flexible functions (e.g. splines on continuous terms, interaction terms)
    # that would reduce model form assumptions (at the cost of some extra variance).
    chisq[i] = g_function_psi(exp.model.formula='radlt0_1 ~ intime + I(intime^2) + cum_rad_lag1 + smoker', t_blip='t_blip', 
                              data=filter(miners, atwork==1), family=binomial())
  }
  tibble(psi=psi.seq, chisq=chisq)
}

gfunction_policy = g_estimation_policy(psi.seq = seq(0, 1, 0.01), occlimit=0.1)
min(gfunction_policy$chisq) # should be close to zero
qplot(x=psi, y=chisq, data=gfunction_policy, geom = 'step', xlab = expression(psi), ylab = expression(chi^2))

(g_estimate_policy <- with(gfunction_policy, psi[which.min(chisq)])) 
lci_policy = with(gfunction_policy, psi[psi<g_estimate_policy][which.min((chisq[psi<g_estimate_policy]-qchisq(.95, 1))^2)])
uci_policy = with(gfunction_policy, psi[psi>g_estimate_policy][which.min((chisq[psi>g_estimate_policy]-qchisq(.95, 1))^2)])

# the 'time ratio' and associated confidence intervals
cat(paste0(signif(exp(g_estimate_policy), 2), " (", signif(exp(lci_policy), 2), ', ', signif(exp(uci_policy), 2), ")"))
# 1.4 (1.1, 1.7)

# Further considerations:
# Administrative censoring (reducing harmful exposures could make some individuals
#   live past the adminstrative censoring date, when in fact, they experienced the event.
#   Artificial censoring is needed in this case, with a re-calculation of the adminstrative
#   censoring date (in the 'unexposed' world) becoming necessary.
# Censoring by loss-to-follow-up: inverse probability of censoring weights (IPCW) are needed
# Competing risks: we may use IPCW to address competing risks, but this changes the 
#   interpretation of the survival ratio
# Last: Functional form of the potential outcome in the survival model:
#   The choice of how the calculated potential outcome enters into the exposure model
#   matters! This is a parametric assumption. When there is administrative censoring,
#   We can use an indicator of death, rather than the actual survival time, which
#   reduces reliance on this assumption, at the cost of efficiency. Other approaches
#   include using relaxed parametric forms For example, we could re-calculate the g-function
#   using a polynomial function of time, as below:


g_function_psi <- function(exp.model.formula, t_blip, data, ...){
  # Alternative way to test independence between t_blip and exposure
  # using a 2 df test under a parametric model with the calculated potential outcome
  #  entered into the exposure model as a quadratic polynomial.
  #  This relaxes the parametric assumptions from the first approach
  f = as.formula(sub("~", paste("~", t_blip, " + I(", t_blip, "^2)", "+"), x = exp.model.formula))
  mod = glm(f, data=data, ...)
  M = summary(mod)$coefficients[2:3,1]
  V = summary(mod)$cov.scaled[2:3, 2:3]
  chisq <- t(M)%*%solve(V)%*%M # chi-squared test statistic for test that two coefs = 0
  chisq
}

# g-estimation with this new function
gfunction2 = g_estimation(psi.seq = seq(0, 1, .01))
# new g-estimate (recall old estimate was 0.59)
(g_estimate2 <- gfunction2$psi[which.min(gfunction2$chisq)])# [1] 0.67
# simulated value is 0.7 (but no guarantee that is the true value in this finite data set)

# estimate may differ due to relaxed parametric assumption, confidence intervals will be wider
qplot(x=psi, y=chisq, data=gfunction2, geom = 'step', xlab = expression(psi), ylab = expression(chi^2))

# determine the confidence intervals by finding the values of psi corresponding to the 
#  closest chi-squared values to the 95th percentile of a chi-squared distribution
#  with the number of df equal to the number of coefficients tested. This is
#  'inverting the test statistic' to find the confidence set
lci2 = gfunction2$psi[gfunction2$psi<g_estimate2][which.min((gfunction2$chisq[gfunction2$psi<g_estimate2]-qchisq(.95, 2))^2)]
uci2 = gfunction2$psi[gfunction2$psi>g_estimate2][which.min((gfunction2$chisq[gfunction2$psi>g_estimate2]-qchisq(.95, 2))^2)]

# the new 'time ratio' and associated confidence intervals (recall original estimate was 0.55 (0.44, 0.67)
cat(paste0(signif(exp(g_estimate2), 2), " (", signif(exp(lci2), 2), ', ', signif(exp(uci2), 2), ")"))
# 2.0 (1.4, 2.3)
# Note point estimates differ while confidence intervals are the same. 
# Inverting a test statistic will not necessarily give bounds that
#  are symmetrical around the point estimate
