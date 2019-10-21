#Analysis of metrical stress free recall experiment

#full record of models run and model selection process
#last edit 14 December 2019 by Amelia 
#impemented with R version 3.4.1 

library(readr)
library(lme4)
library(dplyr)

##############
#Experiment 1#
##############


###############
#FIXED EFFECTS
#Dependent variable: accuracy (correct/wrong) in df: CORR_WRONG(levels: correct, wrong)
#Predictors (levels)
  #FACTORs
      ##Regularity. was the list as a unit regular, or irregular?  
      ###in df: Regularity( levels:irreg,reg)
      ##Pattern.  Was the underlying pattern in the list formed by SW or WS words? 
      ###in df: Pattern (levels: SW,WS)
      ##Target status. Was the individual word in the specified observation 
      ###in target position? in df: TARG_OTH(levels: other, target)
  
  #CONTINUOUS
  ###Number of times marked as "standing out" in the normalization task.
    #in df: "ODD_CHECKS"

################
#RANDOM EFFECTS#
################
#Subject. in df: Subj
#word. in df: word 

exp1dat<-readr::read_csv("Exp1_final_data.csv")
#interpret factors as factors, treatment coding 
  exp1dat$Subj  <- as.factor(exp1dat$Subj)
  exp1dat$Trial <- as.factor(exp1dat$Trial)
  
#exp1 data is not coded for individual word
#create factor for individual word by combining the position variable
  #and the TWord, which specifies the target word for a given list
exp1dat$word <- as.factor(paste(exp1dat$TWord,exp1dat$TListPos,sep=""))

#center numerical variables 
exp1dat$ODD_CHECKS <- exp1dat$ODD_CHECKS-mean(exp1dat$ODD_CHECKS)
exp1dat$FreqLog    <- exp1dat$FreqLog-mean(exp1dat$FreqLog)

exp1mod1 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1+Regularity*Pattern*TARG_OTH|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge

# in order to simply model and achieve convergence, First reduce random effects
#structure (remove interactions in slopes, then reomve slopes,
#then reduce interactions)

exp1mod2 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1+Regularity*Pattern+TARG_OTH|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge

exp1mod3 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1+Regularity+Pattern+TARG_OTH|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge 

exp1mod4 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1+Regularity+Pattern+|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge 

exp1mod5 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1+Regularity|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge

exp1mod6 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH*ODD_CHECKS+
                          (1|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)

#fails to converge 

exp1mod7 <- lme4::glmer(formula = CORR_WRONG~
                          Regularity*Pattern*TARG_OTH+ODD_CHECKS+
                          (1|Subj)+
                          (1|word),
                        family=binomial, 
                        data=exp1dat)
#fails to converge 


###############################
#most complex converging model#
###############################
exp1mod8<-lme4::glmer(formula = CORR_WRONG~
                                Regularity*TARG_OTH+Pattern+ODD_CHECKS+
                                (1|Subj)+
                                (1|word),
                      family = binomial,
                      data = exp1dat)

saveRDS(exp1mod8, "exp1mod8.rds")


##################
#model comparison#
##################
anova(exp1mod8,originalmod1)
#lower AIC,BIC,higher loglik indicates that new model 
#using word as item intercept is better.  


############################################
#Check for participant performance at floor#
############################################

exp1dat$int_accuracy<-recode(exp1dat$CORR_WRONG, "correct" = 1, "wrong" = 0)

by_trial_by_subj_exp1<- exp1dat %>%
                        dplyr::group_by(Subj,Trial) %>%
                        dplyr::summarise(num_corr=sum(int_accuracy))

sum(by_trial_by_subj_exp1$num_corr==0)

#17 out of a total of 3114 trials done across all subjects had zero correct answers. 


##############
#Experiment 2#
##############

#experiment 2 model is  the same as the experiment 1 model, except  there
#is no predictor "odd checks", because experiment 2 involved spliced words 
#and so all tokens of the target word were identical.

###############
#FIXED EFFECTS

#Dependent variable: accuracy (correct/wrong)
#in df: CORR_WRONG(levels: correct, wrong)

#Predictors (levels)
#FACTOR
#### Regularity. was the list as a unit regular, or irregular? 
   #in dataframe (df): Regularity( levels:irreg,reg)
#### Pattern.  Was the underlying pattern in the list formed by SW or WS words?
   #in df: Pattern (levels: SW,WS)
#### Target status. Was the individual word in the specified observation 
   #in target position? in df: TARG_OTH(levels: other, target)


################
#RANDOM EFFECTS#
################
#Subject. in df: Subj
#word. in df: trial 

exp2dat<-readr::read.csv("Exp2_final_data.csv")

exp2dat$SUBJECT<-as.factor(exp2dat$SUBJECT)

#exp2 data is not coded for individual word
##create factor for individual word by combining the position variable
##and the TWord, which specifies the target word for a given list
exp2dat$WORD <- as.factor(paste(exp2dat$TARGET_WORD,exp2dat$TARGET_POSITION,sep=""))

exp2mod1 <- lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                          (1+REG_IRREG*PATTERN*TARG_OTH|SUBJECT)+
                          (1|WORD),
                        family = binomial,
                        data = exp2dat)
#fails to converge 
exp2mod2 <- lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                          (1+REG_IRREG*PATTERN+TARG_OTH|SUBJECT)+
                          (1|WORD),
                        family = binomial,
                        data = exp2dat)
#fails to converge

exp2mod3 <- lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                          (1+REG_IRREG+PATTERN+TARG_OTH|SUBJECT)+
                          (1|WORD),
                        family = binomial,
                        data = exp2dat)
#fails to converge 

exp2mod4 <- lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                          (1+REG_IRREG+PATTERN|SUBJECT)+
                          (1|WORD),
                        family = binomial,
                        data = exp2dat)
#fails to converge 

exp2mod5 <- lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                          (1+REG_IRREG|SUBJECT)+
                          (1|WORD),
                        family = binomial,
                        data = exp2dat)
#fails to converge 

###############################
#most complex converging model#
###############################

exp2mod6 <-lme4::glmer(formula = CORR_WRONG~REG_IRREG*PATTERN*TARG_OTH+
                            (1|SUBJECT)+(1|WORD),
                       family=binomial, 
                       data=exp2dat)
summary(exp2mod6)
saveRDS(exp2mod6, "exp2mod6.rds")


############################################
#Check for participant performance at floor#
############################################

exp2dat$int_accuracy<-recode(exp2dat$CORR_WRONG, "CORRECT" = 1, "WRONG" = 0)

by_trial_by_subj<- exp2dat %>%
                   group_by(SUBJECT,TRIAL) %>% 
                   summarise(num_corr=sum(int_accuracy))

sum(by_trial_by_subj$num_corr==0)

#40 out of a total of 3118 trials across all subjects had zero correct answers. 


