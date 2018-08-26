/**********************************************************************************************************************
* Author: Alex Keil
* Program: workshop_code.sas
* Date: 20180826
* Project: ISEE 2018 workshop: Causal inference foundations and applications in environmental health sciences
* Tasks:
* Data in: miners, coalplant
* Description: Carry out parametric g-formula analysis and strucgtural nested model analysis of data simulated to
*   emulate an occupational study of uranium miners looking at the effects of occupational exposures to radon
*   and all cause mortality, as well as a population study of air pollutants that arise from coal-fired power
*   plants
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/
*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120 PAGESIZE=80 SKIP = 2 FORMDLIM = '-' NOMPRINT NOCENTER; TITLE;TITLE2;FOOTNOTE;


LIBNAME workshop "F:/Talks/Workshops/2018_ISEE_causal/analyses/data/";
%LET outpath = F:/Talks/Workshops/2018_ISEE_causal/analyses/output/;
ODS HTML PATH="&outpath" BODY='saslst.html';
ODS LISTING GPATH = "&outpath";

/********************************************************************************
* Coal plant example - parametric g-formula with several point exposures
*
********************************************************************************/
DATA coalplant;
  SET workshop.coalplant;
RUN;

PROC GENMOD DATA = coalplant;
  TITLE "Parametric g-formula model for point exposures";
  MODEL mdi = as be cd urbanicity as*urbanicity  be*urbanicity  cd*urbanicity  as*be  as*cd  be*cd;
  STORE obsmodel; *the STORE statement allows us to make predictions from this model in a later step;

* Say we know based on local emissions data and discussions with atmospheric
*  chemists that 91% of As, 96% of Be, and 45% of Cd ambient levels come
*  from a local coal-fired power plant;
DATA ints;
  SET coalplant;
  ARRAY a[3] cd as be;
  ARRAY _a[3] _cd _as _be;
  DO i = 1 TO 3;
    *store original values - they will be modified below, but we need to keep them;
    _a[i] = a[i];
  END;
  * no intervention = natural course;
  int = 'NC     ';
  OUTPUT;
  * decommission the coal fired power plant;
  int = 'No coal';
  cd=_cd*(1-0.5);as=_as*(1-0.96);be=_be*(1-0.91);
  OUTPUT;
  DO i = 1 TO 3;
    *restore original values;
    a[i] = _a[i];
  END;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No Cd';
  cd=0;
  OUTPUT;
  cd = _cd;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No As';
  as=0;
  OUTPUT;
  as = _as;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No Be';
  be=0;
  OUTPUT;

* now make predictions based on these new exposures;
PROC PLM RESTORE=obsmodel;
  SCORE DATA = ints OUT=ints PRED = p_mdi;
RUN;

PROC MEANS DATA = ints MEAN;
  TITLE 'G-formula point estimates for mean MDI expected under each intervention';
  CLASS int;
  VAR p_mdi;
RUN;

* get confidence intervals and point estimates for intervention effects by bootstrapping;
%MACRO boot_coal(iter=);
ODS SELECT NONE; OPTIONS NONOTES;
*get sample size;
PROC SQL NOPRINT; SELECT COUNT(mdi) INTO :sampsize FROM coalplant;
%LET i= 0;
%DO %WHILE(%EVAL(&i<=&iter));
%PUT BOOTSTRAP: &i of &iter;
%IF &I=0 %THEN %DO; 
  *0th iteration gets point estimates;
  DATA mcdata; SET coalplant;
%END;
%ELSE %DO; 
  PROC SURVEYSELECT DATA=coalplant OUT=mcdata OUTHITS METHOD=URS N=&sampsize NOPRINT;RUN;
%END;
PROC GENMOD DATA = mcdata;
  TITLE "Parametric g-formula model for point exposures";
  MODEL mdi = as be cd urbanicity as*urbanicity  be*urbanicity  cd*urbanicity  as*be  as*cd  be*cd;
  STORE obsmodel; *the STORE statement allows us to make predictions from this model in a later step;
DATA ints;
  LENGTH int $32;
  SET mcdata;
  ARRAY a[3] cd as be;
  ARRAY _a[3] _cd _as _be;
  DO i = 1 TO 3;
    *store original values - they will be modified below, but we need to keep them;
    _a[i] = a[i];
  END;
  * no intervention = natural course;
  int = 'NC';
  OUTPUT;
  * decommission the coal fired power plant;
  int = 'No coal';
  cd=_cd*(1-0.5);as=_as*(1-0.96);be=_be*(1-0.91);
  OUTPUT;
  DO i = 1 TO 3;
    *restore original values;
    a[i] = _a[i];
  END;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No Cd';
  cd=0;
  OUTPUT;
  cd = _cd;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No As';
  as=0;
  OUTPUT;
  as = _as;
  *eliminate cadmium without intervening on other pollutants;
  int = 'No Be';
  be=0;
  OUTPUT;
PROC PLM RESTORE=obsmodel;
  SCORE DATA = ints OUT=ints PRED = p_mdi;
RUN;
PROC MEANS DATA = ints NOPRINT;
  CLASS int;
  VAR p_mdi;
  OUTPUT OUT=pred_means(KEEP=int p_mdi WHERE=(int ^='')) MEAN=;
RUN;
* get effect estimates;
DATA p_means (KEEP=int p_mdi);
  LENGTH int $32;
  SET pred_means;
  iter = &i;
  RETAIN ref;
  IF int = 'NC' THEN DO;
    int="a) Mean MDI: Natural course";
	ref=p_mdi; *p_mdi will be the mean MDI for natural course OR the intervention effect for other interventions;
    OUTPUT;
  END;
  ELSE IF int = 'No coal' THEN DO;int="b) Effect: Shutdown"; p_mdi=p_mdi-ref; OUTPUT;END;
  ELSE IF int = 'No Cd' THEN DO;int="c) Effect: Attr. MDI Cd"; p_mdi=p_mdi-ref; OUTPUT;END;
  ELSE IF int = 'No As' THEN DO;int="d) Effect: Attr. MDI As"; p_mdi=p_mdi-ref; OUTPUT;END;
  ELSE IF int = 'No Be' THEN DO;int="e) Effect: Attr. MDI Be"; p_mdi=p_mdi-ref; OUTPUT;END;
%IF &I=0 %THEN %DO; 
  *0th iteration contains point estimates from original data;
  DATA pointestimates; SET p_means; RUN; 
%END;
%IF &i=1 %THEN %DO; 
  *be sure there are no results until the 1st iteration;
  DATA bootresults;i=1; PROC SQL NOPRINT; DROP TABLE bootresults; QUIT; 
%END;
PROC APPEND DATA = p_means BASE=bootresults; RUN;
%LET i = %EVAL(&i+1);
%END;
PROC MEANS DATA = bootresults NOPRINT;
 CLASS int;
 VAR p_mdi;
 OUTPUT OUT=boot_se(KEEP=int p_mdi_SE WHERE=(int ^='')) STD=p_mdi_SE;
PROC SORT DATA = pointestimates; BY int;
PROC SORT DATA = boot_se; BY int;
DATA bootsummary;
 MERGE pointestimates boot_se;
 BY int;
 lci = p_mdi-1.96*p_mdi_se;
 uci = p_mdi+1.96*p_mdi_se;
RUN;
TITLE;
ODS SELECT ALL;
  DATA _NULL_;
    SET bootsummary;
    FILE PRINT;
	IF _N_=1 THEN DO;
        PUT "----------------------------------------------------------------------------";
        PUT "| G-formula: policy estimates. Mean MDI, MDI difference                    |";
        PUT "| (compared to natural course) and bootstrap confidence intervals          |";
        PUT "----------------------------------------------------------------------------";
        PUT "|  Intervention/Effect        |    Mean/mean difference (95% CI)           |";
        PUT "----------------------------------------------------------------------------";
	END;
    PUT  @3 int @31 "|" @37 p_mdi 5.3 " (" lci 5.3 ", " uci 5.3 ")" ;
	IF int="a) Mean MDI: Natural course" THEN PUT @31 "|";
  RUN;
  OPTIONS NOTES;
%MEND;
%boot_coal(iter=100);
/********************************************************************************
* Uranium miners examples
*
********************************************************************************/


DATA miners;
  SET workshop.miners;
  *some variables for convenience later;
  IF atwork=1 THEN lograd = LOG(rad);
RUN;



/********************************
* estimation of policy effects with the parametric g-formula
********************************/

  /* modeling: we do this in a macro because later it will be easier for bootstrapping*/
%MACRO get_coefs(indata);
  PROC GENMOD DATA = &indata;
    TITLE 'Exposure model for g-formula estimates of occupational-radon policy effects on mortality';
    WHERE atwork=1;
    MODEL lograd = outtime cum_rad_lag1 smoker / D=NORMAL LINK=ID CONVERGE=1E-16;
	ODS OUTPUT ParameterEstimates = ecoefs (KEEP=parameter estimate);
	*because this is a linear model, the last coefficient will be the scale (~ sd of the residuals);
  PROC GENMOD DATA = &indata DESCENDING;
    TITLE 'Employment (confounder) model for g-formula estimates of occupational-radon policy effects on mortality';
    WHERE atwork=1 or leavework=1;
    MODEL leavework = outtime outtime*outtime rad_lag1 cum_rad_lag1 smoker  / D=B LINK=LOGIT;
	ODS OUTPUT ParameterEstimates = wcoefs (KEEP=parameter estimate WHERE=(parameter ^= "Scale"));
  PROC GENMOD DATA = &indata DESCENDING;
    TITLE 'Death (outcome) model for g-formula estimates of occupational-radon policy effects on mortality';
    MODEL dead = outtime atwork cum_rad cum_rad*cum_rad smoker smoker*cum_rad / D=B LINK=LOGIT CONVERGE=1E-16;
	ODS OUTPUT ParameterEstimates = dcoefs (KEEP=parameter estimate WHERE=(parameter ^= "Scale"));
  RUN;
  /* global macro variables that contain the coefficients and counts of coefficients */
  %GLOBAL ecoefs necoefs wcoefs nwcoefs dcoefs ndcoefs;
  PROC SQL NOPRINT;
    SELECT COUNT(estimate), estimate FORMAT=best32. INTO :necoefs, :ecoefs SEPARATED BY " "  FROM ecoefs;
    SELECT COUNT(estimate), estimate FORMAT=best32. INTO :nwcoefs, :wcoefs SEPARATED BY " "  FROM wcoefs;
    SELECT COUNT(estimate), estimate FORMAT=best32. INTO :ndcoefs, :dcoefs SEPARATED BY " "  FROM dcoefs;
  QUIT;
  TITLE;
%MEND;

%MACRO gformula_risks(indata, outdata, endtime=20, mc_iter=50000, seed=12312);
PROC SURVEYSELECT DATA=miners(WHERE=(intime=0)) 
   OUT=pseudo_cohort_bl OUTHITS N=&mc_iter METHOD=URS NOPRINT;

DATA &outdata(DROP=_E: _W: _D: mul mue mud lograd);
 SET pseudo_cohort_bl;
 * bring in model coefficients using SAS arrays;
 ARRAY _e[&necoefs] (&ecoefs);
 ARRAY _w[&nwcoefs] (&wcoefs);
 ARRAY _d[&ndcoefs] (&dcoefs);
 gfid = _N_;
 CALL STREAMINIT(&seed);
 DO intervention = -1, 0, 1, 3; *4 interventions: natural course, unexposed, limit at 0.1, limit at 0.3;
   *initial values for time-varying variables;
   atwork=1;
   durwork=0;
   durwork_lag1 = 0;
   leavework=0;
   rad_lag1 = 0;
   cum_rad_lag1 = 0;
   dead = 0;
   cum_rad = 0;
   DO intime = 0 TO (&ENDTIME-1);
     outtime=intime+1;
	 IF atwork THEN DO;
	   IF intime>0 THEN DO;
	     *in observed data, no one leaves in first year - enforce that here;
	     *log-odds of leaving work;
         mul = _w[1] + _w[2]*outtime + _w[3]*outtime*outtime + _w[4]*rad_lag1 + _w[5]*cum_rad_lag1 + _w[6]*smoker;
         *random draw from a bernoulli distribution with probability of leaving work given 
		  by the inverse logit transform of the log-odds of leaving work;
         leavework = RAND('bernoulli', 1/(1+exp(-mul))); 
       END;
	   IF leavework = 1 THEN atwork = 0;
	   durwork = durwork+atwork;
	 END; *atwork;
	 * exposure;
	 /*********
	 * static intervention: always unexposed
	 **********/
	 IF intervention = 0 THEN rad = 0;
	 ELSE IF intervention IN (-1 1 3) THEN DO;
	   /*********
	   * dynamic/policy interventions: if at work, exposed below some limit (natural course kept below the observed limit);
	   **********/
       IF atwork THEN DO;
	     *mean log exposure;
         mue = _e[1] + _e[2]*outtime  + _e[3]*cum_rad_lag1  + _e[4]*smoker ;
         *log exposure is a random draw from the regression mean and std. dev. of the residuals, capped at empirical maximum;
         DRAWLOGEXPOSURE: lograd = MIN(RAND('NORMAL', mue, _e[5]), LOG(7.321276));
         rad = EXP(lograd);
		 *this antiquated bit of programming uses 'goto' statements to take another draw from the log-exposure distribution if it is above the limit;
		 * this is one way to draw from a truncated normal distribution with an upper bound at the (log) occupational limit under consideration;
		 IF intervention = 1 AND rad>0.1 THEN GOTO DRAWLOGEXPOSURE;
		 ELSE IF intervention = 3 AND rad>0.3 THEN GOTO DRAWLOGEXPOSURE;
	   END; *atwork;
	   ELSE DO; lograd=.; rad = 0; END;
	 END; *interventions;
     cum_rad = cum_rad+rad;
	 /*********
	 * discrete hazard of death
	 **********/
     mud = _d[1] + _d[2]*outtime + _d[3]*atwork + _d[4]*cum_rad + _d[5]*cum_rad*cum_rad + _d[6]*smoker + _d[7]*smoker*cum_rad;
	 pdead = 1/(1+exp(-mud));
	 *update cumulative incidence with product limit estimator;
	 IF intime = 0 THEN cum_incidence = pdead;
     ELSE cum_incidence = cum_incidence + (1-cum_incidence)*pdead;
	 /*
	 *alternative way to get cumulative incidence with the Breslow estimator;
	 IF intime = 0 THEN cum_hazard = pdead;
     ELSE cum_hazard = cum_hazard + pdead;
     cum_incidence = 1-EXP(-cum_hazard);
	 */
     OUTPUT;
	 *leavework should only equal 1 in a single observation;
	 IF leavework=1 THEN leavework = 0;
	 *lagged variables;
	 rad_lag1 = rad;
	 cum_rad_lag1 = cum_rad;
	 durwork_lag1 = durwork;
   END;*DO intime = 0 TO (&ENDTIME-1);
 END;*DO intervention = -1, 0, 1, 3;
RUN;
%MEND;
%MACRO gformula_effectestimates(indata, outdata, endtime=20);
  PROC MEANS DATA = &INDATA NOPRINT;
    CLASS intervention outtime;
	VAR cum_incidence;
	OUTPUT OUT=cilong(WHERE=(intervention>.z AND outtime>.z) DROP=_:) MEAN=;
  RUN;
  DATA &outdata (KEEP=outtime ci: rd:);
   SET cilong END=EOF;
   ARRAY ci[4];
   ARRAY _ci[4, &endtime] _TEMPORARY_;
   ARRAY int[4] (-1 0 1 3);
   DO i = 1 TO 4;
     IF intervention = int[i] THEN _ci[i, outtime] = cum_incidence;
   END;
   IF eof THEN DO;
     DO outtime = 0 TO &endtime;
       IF outtime=0 THEN DO M = 1 TO 4;
         ci[m]=0;;
       END;
       ELSE DO M = 1 TO 4;
         ci[m]=_ci[m, outtime];;
       END;
	   rd_unex = ci[2]-ci[1];
	   rd_exp1 = ci[3]-ci[1];
	   rd_exp3 = ci[4]-ci[1];
	   RENAME ci1 = ci_nc ci2=ci_unex ci3=ci_exp1 ci4=ci_exp3;
       OUTPUT;
	 END;
   END;
 RUN;
%MEND;

%get_coefs(miners);
*PROC PRINT DATA = ecoefs NOOBS;
*  TITLE "example of dataset holding coefficients";
*RUN;
%gformula_risks(indata=miners, outdata=pseudo_cohort, mc_iter=100000, endtime=20);
%gformula_effectestimates(indata=pseudo_cohort, outdata=pointestimates, endtime=20);


* time ratio for comparison with SNM;
PROC SQL;
  TITLE 'Policy time ratio for comparison with structural nested model (SNM TR=1.36)';
  CREATE TABLE a AS 
     SELECT
     SUM(1-cum_incidence) AS elifnc FROM pseudo_cohort(WHERE=(intervention=-1));
  CREATE TABLE b AS 
     SELECT
     SUM(1-cum_incidence) AS elifexp1 FROM pseudo_cohort(WHERE=(intervention=1));
  SELECT elifexp1/elifnc AS TR FROM a, b;
QUIT;
* 1.30;

*getting observed cumulative incidence;
PROC PHREG DATA = miners NOPRINT;
  MODEL (intime outtime)*dead(0)=;
  BASELINE OUT=obs_surv SURVIVAL=s_observed;
RUN;

DATA gfplotdata;
  SET pointestimates obs_surv;
  ci_observed = 1-s_observed;
  LABEL ci_observed = "Observed";
  LABEL ci_nc = "GF: Natural course";
  LABEL ci_unex = "GF: Never exposed";
  LABEL ci_exp1 = "GF: Limit 0.1";
  LABEL ci_exp3 = "GF: Limit 0.3";

ODS LISTING GPATH = "&outpath";
ODS GRAPHICS ON / RESET = ALL IMAGENAME="sas_gformulaCI";
ODS ESCAPECHAR='~';
PROC SGPLOT DATA = gfplotdata;
  TITLE 'G-formula estimates of cumulative incidence under alternative polcies on annual radon exposure';
  SERIES x = outtime y = ci_observed;
  SERIES x = outtime y = ci_nc;
  SERIES x = outtime y = ci_unex;
  SERIES x = outtime y = ci_exp1;
  SERIES x = outtime y = ci_exp3;
  YAXIS LABEL = "Cumulative incidence";
  XAXIS LABEL = "Time";
QUIT;
TITLE;
* getting confidence intervals for risk difference at time 10;
%MACRO BOOTSTRAP_GF(indata, iter=10, endtime=10);
  ODS EXCLUDE ALL; ODS NORESULTS;
  OPTIONS NONOTES;
  DATA boot_samples; a=1;
  PROC SQL NOPRINT;
    SELECT COUNT (UNIQUE(ID)) INTO :sampsize FROM &indata;
    DROP TABLE boot_samples;
  %LET i = 1;
  %DO %WHILE(%EVAL(&i<=&iter));
  %PUT BOOTSTRAP_GF: iteration &i of &iter;
    PROC SURVEYSELECT DATA=&indata OUT=bootpop METHOD=URS N=&sampsize NOPRINT; SAMPLINGUNIT id;RUN;
    %get_coefs(bootpop);
    %gformula_risks(indata=bootpop, outdata=pseudo_cohort, endtime=&endtime, mc_iter=50000, seed=&i);
    %gformula_effectestimates(indata=pseudo_cohort, outdata=bootpe, endtime=&endtime);
	PROC APPEND DATA=bootpe(WHERE=(outtime=&endtime)) BASE=boot_samples;
  %LET i = %EVAL(&I+1);
  %END;

  PROC MEANS DATA = boot_samples NOPRINT;
    VAR ci_: rd:;
	OUTPUT OUT=bootsum STD=/AUTONAME;
  DATA gf_policy;
    MERGE pointestimates(WHERE=(outtime=&endtime)) bootsum;
	  int = "Limit = 0.3    ";
	  rd = rd_exp3;
	  lci = rd-1.96*rd_exp3_stddev;
	  uci = rd+1.96*rd_exp3_stddev;
	  OUTPUT;
	  int = "Limit = 0.1    ";
	  rd = rd_exp1;
	  lci = rd-1.96*rd_exp1_stddev;
	  uci = rd+1.96*rd_exp1_stddev;
	  OUTPUT;
	  int = "Never exposed";
	  rd = rd_unex;
	  lci = rd-1.96*rd_unex_stddev;
	  uci = rd+1.96*rd_unex_stddev;
	  OUTPUT;
  RUN;
  ODS EXCLUDE NONE; ODS RESULTS;
  OPTIONS NOTES;

  DATA _NULL_;
    SET gf_policy END=eof;
    FILE PRINT;
	IF _N_=1 THEN DO;
        PUT "----------------------------------------------------------------------------";
        PUT "| G-formula: policy estimates. &endtime-year Risk difference estimates           |";
        PUT "| (compared to natural course) and bootstrap confidence intervals          |";
        PUT "----------------------------------------------------------------------------";
        PUT "|  Intervention   |    Risk difference (95% CI)                            |";
        PUT "----------------------------------------------------------------------------";
	END;
        PUT  @3 int @19 "|" @25 rd 5.3 " (" lci 5.3 "," uci 5.3 ")" ;
  RUN;
%MEND;

*this takes a while;
*%BOOTSTRAP_GF(indata=miners, iter=200, endtime=10);
*just do a few iterations to get the idea;
%BOOTSTRAP_GF(indata=miners, iter=10, endtime=10);


/********************************************************************************
* structural nested accelerated failure time model (fit with g-estimation)
* an SNM is a model for the effect of a 'blip' of treatment, conditional on being untreated in the future
* it compares the outcomes we would observe under no treatment to the outcomes we would observe
*   if treated now and then kept untreated in the future. This is an admittedly odd
*   way to think about causal effects, but it corresponds to a couple things we might be familiar with:
*   1) in time fixed settings, the only 'blip' is the time-fixed exposure so a standard (say) linear
*      model for a continuous outcome, adjusted for covariates is equivalent to a conditional 
*      structural nested model
*   2) in time-varying exposure settings, if we add all of the blip effects for "always treated"
*      we estimate the same "always treated" counterfactual as we do in other causal methods.
*      In other words, the "blip effect" is just a way of decomposing the "always" treated effect
*      
********************************************************************************/




* setting up the data;
DATA gestimation;
  SET miners;
  DO psi = 0 TO 1 BY 0.01;
    *blipping down;
  * blip function: the effect of a blip of exposure, holding future exposure
  *  to equal zero. Calculated for each row of the dataset.
  * log-linear SNM
  * observed "duration" of time for observation (e.g. 1 person year)
  * is divided by exp(X*psi), which gives the expected 
  *  contribution to the total survival time;
    pyr0 = EXP(psi*cum_rad)*(outtime-intime);
    OUTPUT;
  END;


PROC SORT DATA = gestimation;
  BY psi id DESCENDING intime;

*calculationg
  * calculate total expected lifetime under no exposure from time 'start'
  * until the end of the observed time. Calculated for each observation
  * in the dataset
  * the method described in the theory:;
DATA gestimation;
  SET gestimation;
  BY psi id DESCENDING intime;
  RETAIN revcumpyr0;
  IF first.id THEN revcumpyr0=0;
  revcumpyr0 = revcumpyr0 + pyr0;
  t_blip = revcumpyr0 + intime;
RUN;

  * a slightly less efficient method that is easier in most software
  * this calculates the potential outcome under "never exposed" for all observations
  * See Hernan et al (2005) Pharmacoepidem Dr S:;
PROC SQL;
  CREATE TABLE gestimation2 AS 
    SELECT *, SUM(pyr0) AS t_blip2 FROM gestimation GROUP BY psi, id;
QUIT;


  /* calculate the 'g-function', or the test-statistic comparing observed
  *  exposure with the blipped-down time (this is a function of a particular
  *  value of psi and does not estimate psi!);
    * note that the exposure model is fit only to the person time for which individuals
    *  are at work! This leverages the "selective ignorability" condition, where
    *  we need only test exchangeability in a subset of the data in which it is 
    *  expected to hold. In this case, we necessarily fit the model to the employed
    *  person time because exposure does not occur outside of work (i.e. non-positivity)
    * Model: pooled linear model for log(exposure) as a function of time (quadratic), 
    *  cumulative prior exposure, employment duration, and smoking. In practice, we would usually fit
    *  more complex, flexible functions (e.g. splines on continuous terms, interaction terms)
    * that would reduce model form assumptions (at the cost of some extra variance).*/

ODS SELECT NONE;
PROC GENMOD DATA = gestimation2;
  TITLE 'Exposure model for g-estimation of log-time ratio per unit of cumulative radon exposure';
  BY psi;
  WHERE atwork=1;
  MODEL lograd = t_blip intime intime*intime cum_rad_lag1 smoker;
  ODS OUTPUT ParameterEstimates = gfunction(KEEP=chisq psi parameter WHERE=(parameter="t_blip"));
  *below is code that uses the less-efficient estimator;
  *MODEL lograd = t_blip2 intime intime*intime cum_rad_lag1 smoker;
  *ODS OUTPUT ParameterEstimates = gfunction(KEEP=chisq psi parameter WHERE=(parameter="t_blip2"));
RUN;
ODS SELECT ALL;

PROC SQL NOPRINT;
  CREATE TABLE gfunction2 AS SELECT psi, chisq, MIN(chisq) AS minchi FROM gfunction;

/* The g-function is a series of 'objective function' (chi-squared values) 
*  across a grid of possible values of psi - this is a 
*  quick and dirty way to get g-estimates. We could use
*  an optimization routine (e.g. nelder-mead simplex algorithm)
*  to find the minimum, but a grid search over a set of plausible values is a 
*  more sure way to find the minimum since the objective function may have 
*  local minima - this is challenging when there are multiple psi parameters,
*  so optimization may be a better option in that case
* Here: we get chi-squared values for a grid of possible values of psi ranging
*  from -1 to -0.4 (the data are simulated, so I know the true estimate lies in this
*  range. In practice, we might want to have a larger range.*/

/* Plotting the g-function. the g-estimate of psi is given by the value of psi for
* which the following hold:
*  1) the calculated potential outcomes (given by the structural model) are 
*    statistically independent from the observed exposure (given the exposure model)
*  2) the p-value for the association between the calculated potential outcomes and the
*     observed exposure is 1.0, in a test of association that is adjusted for
*     time-varying and baseline confounders
*  or 3) the coefficient for a variable with the calculated potential outcomes is zero
*     in a model for the exposure, controlling for time-varying and baseline confounders
*  These conditions are stated in order of specificity: we could establish independence
*   using (for example) log-rank tests. In this example, we use a regression model for the
*   log(exposure) [NOT cumulative exposure], adjusted for time-fixed and time-varying confounders.
*   Note that this is simply an easy (relatively) to understand approach and is not 
*   the most computationally efficient!
*   In the exposure model, we use a Wald statisic (coefficient^2/variance) ~ Chi^2(df=1) to test
*   whether the coefficient for the (guessed) potential outcomes are independent of observed
*   exposure. By conditional exchangeability, the calculated potential outcomes and the exposure
*   will be independent at the true value of psi. A confidence interval set is created by
*   inverting the Chi-squared statistic to find the points closest to the g-estimate of psi
*   at which the Chi-squared statistic is equal to the critical value of 3.84 (we could be
*   more exact and use linear interpolation between the grid-points, but an approximation will 
*   suffice here).*/

DATA psiest(KEEP=psihat lci uci tr:);
  SET gfunction2 END=EOF;
  RETAIN reachedmin psihat uci lci 0 chidistu chidistl 1e6;
  IF chisq = minchi THEN DO;
    psihat = psi;
	CALL SYMPUT('psihat', PUT(psi, BEST9.));
	reachedmin = 1;
  END;
  IF reachedmin=0 AND (chisq-3.84)**2<=chidistl THEN DO;
    lci = psi;
    chidistl=(chisq-3.84)**2;
  END;
  IF reachedmin=1 AND (chisq-3.84)**2 < chidistu THEN DO;
    uci = psi;
    chidistu=(chisq-3.84)**2;
  END;
  TR = EXP(psihat);
  TRlci = EXP(lci);
  TRuci = EXP(uci);
  IF eof THEN OUTPUT;

ODS GRAPHICS ON / RESET = ALL IMAGENAME="sas_gfunction";
ODS ESCAPECHAR='~';
PROC SGPLOT DATA = gfunction;
  TITLE 'G-estimate of psi and 95% confidence intervals';
  STEP x = psi y = chisq;
  REFLINE 3.84 / AXIS = y LABEL="95% CI";
  REFLINE &psihat / AXIS = x LABEL="~{unicode psi}~{unicode hat}";
  YAXIS LABEL = '~{unicode chi}^2';
  XAXIS LABEL = "~{unicode psi}";
QUIT;
TITLE;
  DATA _NULL_;
    SET psiest END=eof;
    FILE PRINT;
        PUT "----------------------------------------------------------------------------";
        PUT "| G-estimation: structural nested accelerated failure time models          |";
        PUT "| psi (log time ratio) and time ratio estimates and confidence intervals   |";
        PUT "----------------------------------------------------------------------------";
        PUT "|         Psi (95% CI)                |           Time ratio (95% CI)      |";
        PUT "----------------------------------------------------------------------------";
        PUT @10 psihat 4.2 " (" lci 4.2 "," uci 4.2 ")" @50 TR 4.2 " (" trlci 4.2 "," truci 4.2 ")";
  RUN;

/********************************
* g-estimation of policy effects
********************************/


  * setting up the data;
DATA gestimation_policy;
  SET miners;
  DO psi = 0 TO 1.0 BY 0.01;
  * the 'policy exposure' - you are compliant with the policy (referent)
   if you are below the occupational limit of 0.1, OR if you are not at work,
  Thus, the occupational policy we are testing is "If at work, remain unexposed";
  radlt0_1 = (rad>0.1);

   *blipping down;
    * note that here we define the causal effect of interest;
  * log-linear SNM;
    pyr0 = EXP(psi*radlt0_1)*(outtime-intime);
    OUTPUT;
  END;


PROC SORT DATA = gestimation_policy;
  BY psi id DESCENDING intime;

*calculationg
  * calculate total expected lifetime under no exposure from time 'start'
  * until the end of the observed time. Calculated for each observation
  * in the dataset
  * the method described in the theory:;
DATA gestimation_policy2;
  SET gestimation_policy;
  BY psi id DESCENDING intime;
  RETAIN revcumpyr0;
  IF first.id THEN revcumpyr0=0;
  revcumpyr0 = revcumpyr0 + pyr0;
  t_blip = revcumpyr0 + intime;
RUN;
ODS SELECT NONE;
PROC GENMOD DATA = gestimation_policy2;
  TITLE 'Exposure model for g-estimation of log-time ratio for policy to limit annual radon to <0.1 units';
  BY psi;
  WHERE atwork=1;
  MODEL radlt0_1 = t_blip intime intime*intime cum_rad_lag1 smoker / D=B LINK=LOGIT;
  ODS OUTPUT ParameterEstimates = gfunction_policy(KEEP=chisq psi parameter WHERE=(parameter="t_blip"));
RUN;
TITLE;
ODS SELECT ALL;

PROC SQL NOPRINT;
  CREATE TABLE gfunction_policy2 AS SELECT psi, chisq, MIN(chisq) AS minchi FROM gfunction_policy;

DATA psiest_policy(KEEP=psihat lci uci tr:);
  SET gfunction_policy2 END=EOF;
  RETAIN reachedmin psihat uci lci 0 chidistu chidistl 1e6;
  IF chisq = minchi THEN DO;
    psihat = psi;
	CALL SYMPUT('psihatp', PUT(psi, BEST9.));
	reachedmin = 1;
  END;
  IF reachedmin=0 AND (chisq-3.84)**2<=chidistl THEN DO;
    lci = psi;
    chidistl=(chisq-3.84)**2;
  END;
  IF reachedmin=1 AND (chisq-3.84)**2 < chidistu THEN DO;
    uci = psi;
    chidistu=(chisq-3.84)**2;
  END;
  TR = EXP(psihat);
  TRlci = EXP(lci);
  TRuci = EXP(uci);
  IF eof THEN OUTPUT;

  DATA _NULL_;
    SET psiest_policy END=eof;
    FILE PRINT;
        PUT "----------------------------------------------------------------------------";
        PUT "| G-estimation: policy structural nested accelerated failure time models   |";
        PUT "| psi (log time ratio) and time ratio estimates and confidence intervals   |";
        PUT "----------------------------------------------------------------------------";
        PUT "|         Psi (95% CI)                |           Time ratio (95% CI)      |";
        PUT "----------------------------------------------------------------------------";
        PUT @10 psihat 4.2 " (" lci 4.2 "," uci 4.2 ")" @50 TR 4.2 " (" trlci 4.2 "," truci 4.2 ")";
  RUN;

/* Further considerations:
* Administrative censoring (reducing harmful exposures could make some individuals
*   live past the adminstrative censoring date, when in fact, they experienced the event.
*   Artificial censoring is needed in this case, with a re-calculation of the adminstrative
*   censoring date (in the 'unexposed' world) becoming necessary.
* Censoring by loss-to-follow-up: inverse probability of censoring weights (IPCW) are needed
* Competing risks: we may use IPCW to address competing risks, but this changes the 
*   interpretation of the survival ratio
* Last: Functional form of the potential outcome in the survival model:
*   The choice of how the calculated potential outcome enters into the exposure model
*   matters! This is a parametric assumption. When there is administrative censoring,
*   We can use an indicator of death, rather than the actual survival time, which
*   reduces reliance on this assumption, at the cost of efficiency. Other approaches
*   include using relaxed parametric forms For example, we could re-calculate the g-function
*   using a polynomial function of time, as below:*/


ODS SELECT NONE;
PROC GENMOD DATA = gestimation2;
  TITLE 'Alternative exposure model for g-estimation of log-time ratio per unit of cumulative radon exposure';
  BY psi;
  WHERE atwork=1;
  MODEL lograd = t_blip t_blip*t_blip intime intime*intime cum_rad_lag1 smoker / COVB;
  ODS OUTPUT ParameterEstimates = gfunctionb(KEEP=psi estimate parameter WHERE=(parameter="t_blip" OR parameter="t_blip*t_blip"))
              COVB=covmat(KEEP=prm2 prm3 rowname psi WHERE=(rowname='Prm2' OR rowname='Prm3'));
RUN;
ODS SELECT ALL;

PROC IML;
  TITLE 'Calculation of chi-squared statistic for g-estimation under alternative exposure model';
  USE gfunctionb();
  READ all var{psi estimate} INTO beta;
  USE covmat();
  READ all var{psi prm2 prm3} INTO cov;
  CLOSE gfunction;
  unique_rows = UNIQUEBY(beta, 1,1:NROW(beta));
  chisq = J(NROW(unique_rows), 1);
  psi = J(NROW(unique_rows), 1);
  DO p=1 TO NROW(unique_rows);
	IF p=nrow(unique_rows) then index=unique_rows[p]:nrow(beta);
	ELSE index = unique_rows[p]:unique_rows[p+1]-1;
	B = beta[index,2];
	V = cov[index,2:3];
	psi[p] = beta[unique_rows[p],1];
	chit = T(B)*INV(V)*B;
	chisq[p] = T(B)*INV(V)*B;
  END;
  CREATE gfunctionb2  VAR{psi chisq};
  APPEND;
  CLOSE gfunctionb2;
QUIT;
PROC SQL NOPRINT;;
  CREATE TABLE gfunctionb3 AS SELECT psi, chisq, MIN(chisq) AS minchi FROM gfunctionb2;

DATA psiestb(KEEP=psihat lci uci tr:);
  SET gfunctionb3 END=EOF;
  RETAIN reachedmin psihat uci lci 0 chidistu chidistl 1e6;
  IF chisq = minchi THEN DO;
    psihat = psi;
	CALL SYMPUT('psihat', PUT(psi, BEST9.));
	reachedmin = 1;
  END;
  IF reachedmin=0 AND (chisq-3.84)**2<=chidistl THEN DO;
    lci = psi;
    chidistl=(chisq-3.84)**2;
  END;
  IF reachedmin=1 AND (chisq-3.84)**2 < chidistu THEN DO;
    uci = psi;
    chidistu=(chisq-3.84)**2;
  END;
  TR = EXP(psihat);
  TRlci = EXP(lci);
  TRuci = EXP(uci);
  IF eof THEN OUTPUT;


  DATA _NULL_;
    SET psiestb END=eof;
    FILE PRINT;
        PUT "----------------------------------------------------------------------------";
        PUT "| G-estimation: structural nested accelerated failure time models          |";
        PUT "| psi (log time ratio) and time ratio estimates and confidence intervals   |";
        PUT "----------------------------------------------------------------------------";
        PUT "|         Psi (95% CI)                |           Time ratio (95% CI)      |";
        PUT "----------------------------------------------------------------------------";
        PUT @10 psihat 4.2 " (" lci 4.2 "," uci 4.2 ")" @50 TR 4.2 " (" trlci 4.2 "," truci 4.2 ")";
  RUN;




ODS HTML close;
