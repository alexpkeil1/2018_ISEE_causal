**** PROPENSITY SCORES ********

* Thanks to Sam Harper (2013) for code used in McGill 672 class

use "C:\Users\jkaufm5\Dropbox\JSK\Invited Lectures\Michigan Summer 2018\Wednesday\nswre74.dta"

* label some variables
label var treat "treatment group"
label var age "age (years)"
label var ed "years of education"
label var black "black race"
label var hisp "hispanic ethnicity"
label var married "married"
label var nodeg "no high school degree"
label var re74 "1974 earnings"
label var re75 "1975 earnings"
label var re78 "1978 earnings"
label var age2 "age squared"
label define noyes 0 no 1 yes, modify
label values treat black hisp married nodeg noyes

* means of pre-treatment covariates 

* table of means of covariates across treated/control groups
tabstat age age2 ed black hisp nodeg married re74 re75, by(treat) stat(mean sd) format(%8.0g)

* balance by hand
foreach var of varlist age age2 ed black hisp nodeg married re74 re75 {
      qui sum `var' if treat==1
      scalar m`var'_t = r(mean)
      scalar sd`var'_t = r(sd)
      disp as text "Mean(SD) of `var' for treated: " as result %6.2f m`var'_t " (" %6.2f sd`var'_t ")"
      qui sum `var' if treat==0
      scalar m`var'_c = r(mean)
      scalar sd`var'_c = r(sd)
      disp as text "Mean(SD) of `var' for control: " as result %6.2f m`var'_c " (" %6.2f sd`var'_c ")"
      scalar d`var' = (100*abs(m`var'_t - m`var'_c)) / sqrt(sd`var'_t^2 )
      disp as text "Standardized difference for `var' = " as result %3.1f d`var'
      disp " "
      }
 
* build propensity score models for nsw ("by hand")

* predict treatment status
logit treat age age2 ed black hisp nodeg married re75 re74

* output predicted probability of treatment
predict ps

* summary statistics for propensity score, by treatment
bysort treat: sum ps, det
histogram ps, by(treat)

* nearest neighbor matching (1:1) without replacement

* randomly order data in case match ties
set seed 123456 // set seed to reproduce results

gen ranorder = runiform() 
sort ranorder

* create propensity score
psmatch2 treat age age2 ed black hisp nodeg married re75 re74, logit neighbor(1) noreplacement

* evaluate balance (Rosenbaum and Rubin formula)
pstest age age2 ed black hisp nodeg married re75 re74, both graph

* checking balance by hand (Stuart 2010 formula)
foreach var of varlist age age2 ed black hisp nodeg married re75 re74 {

	qui sum `var' if _weight==1 & _treat==1
	scalar m`var'_t = r(mean)
	scalar sd`var'_t = r(sd)
	disp as text "Mean(SD) of `var' for treated: " as result %6.2f m`var'_t " (" %6.2f sd`var'_t ")"
	
	qui sum `var' if _weight==1 & _treat==0
	scalar m`var'_c = r(mean)
	scalar sd`var'_c = r(sd)
	disp as text "Mean(SD) of `var' for control: " as result %6.2f m`var'_c " (" %6.2f sd`var'_c ")"
	
	scalar d`var' = (100*abs(m`var'_t - m`var'_c)) / sqrt((sd`var'_t^2))
	
	disp as text "Standardized difference for `var' = " as result %3.1f d`var'
	disp " "
	
	}
	 
* modified nearest neighbor matching (1:1) without replacement (optional)

* generate squared term for 1974 earnings
gen re74_2 = re74*re74
label var re74_2 "1974 earnings squared"

* randomly order data in case match ties
replace ranorder = runiform() 
sort ranorder

* create propensity score
psmatch2 treat age age2 ed black hisp nodeg married re75 re74 re74_2, logit neighbor(1) noreplacement

* evaluate balance 
pstest age age2 ed black hisp nodeg married re75 re74, both graph  

* checking balance by hand 
foreach var of varlist age age2 ed black hisp nodeg married re75 re74 {

	qui sum `var' if _weight==1 & _treat==1
	scalar m`var'_t = r(mean)
	scalar sd`var'_t = r(sd)
	disp as text "Mean(SD) of `var' for treated: " as result %6.2f m`var'_t " (" %6.2f sd`var'_t ")"
	
	qui sum `var' if _weight==1 & _treat==0
	scalar m`var'_c = r(mean)
	scalar sd`var'_c = r(sd)
	disp as text "Mean(SD) of `var' for control: " as result %6.2f m`var'_c " (" %6.2f sd`var'_c ")"
	
	scalar d`var' = (100*abs(m`var'_t - m`var'_c)) / sqrt((sd`var'_t^2))
	
	disp as text "Standardized difference for `var' = " as result %3.1f d`var'
	disp " "
	
	}

* treatment effect in final matched sample

regress re78 treat if _weight==1
regress re78 treat age age2 ed black hisp nodeg married re75 re74 if _weight==1

* one-step propensity matched treatment effect
psmatch2 treat age age2 ed black hisp nodeg married re75 re74 re74_2, logit neighbor(1) noreplacement outcome(re78)

* examine support
tab _treat _support
psgraph, name(pscore1, replace)    /* Examine support graphically */

* bootstrapped standard error for ATT:
bootstrap r(att), reps(500): psmatch2 treat age age2 ed black hisp nodeg married re75 re74 re74_2, logit neighbor(1) noreplacement outcome(re78)

* using Stata treatment effects suite of commands

teffects psmatch (re78) (treat age ed black hisp married nodeg re74 age2 re74_2), atet
teffects psmatch (re78) (treat age ed black hisp married nodeg re74 age2 re74_2)
teffects psmatch (re78) (treat age ed black hisp married nodeg re74 age2 re74_2), nneighbor(2) nopvalues cformat(%9.2f) pformat(%5.2f) sformat(%8.2f) caliper(.2)

* post-estimation commands
teffects overlap
tebalance summarize, baseline
tebalance summarize
tebalance density ed

* easier to specify more flexible model 
teffects psmatch (re78) (treat age ed black hisp married nodeg c.age#(c.age c.ed i.nodeg) c.re74#(c.re74 i.black))
tebalance summarize
tebalance box ed







	
	
	
	
	
	



















