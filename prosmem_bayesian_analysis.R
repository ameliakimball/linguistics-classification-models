
rm(list=ls())

# loading the data frame with the results
data <- read.csv("Final_data_Aug_15_step_pairs.csv")
summary(data)

# make subject a factor
data$SUBJ_SHORT <- as.factor(data$SUBJ_SHORT)
data$TIME_LAG<-as.factor(data$TIME_LAG)
data$STEP_PAIR_COLL<-as.factor(data$STEP_PAIR_COLL)

# transforming the dv into a numeric (1/0)
data$dv <- ifelse(data$CORR_RESP=="CORR",1, ifelse(data$accuracy=="WRONG",0,NA))

#set up forward difference coding for time lag 
my.forward.diff.4 = matrix(c(3/4, -1/4, -1/4, -1/4, 1/2, 1/2, -1/2, -1/2,
                             1/4, 1/4, 1/4, -3/4), ncol = 3)
#apply forward difference coding to the contrasts
contrasts(data$TIME_LAG)<-my.forward.diff.4



# create a df with only different data
diff_data<-subset(data,data$CORR_RESP=="DIFF")

#interpret timelag as a factor
diff_data$TIME_LAG<-as.factor(diff_data$TIME_LAG)
table(diff_data$TIME_LAG,diff_data$USER_CORR)


#drop unused levels in all variables to be used in model ( STEP_PAIR_COLL, SUBJ_SHORT,TIME_LAG, EXPERIMENT)
diff_data$CORR_RESP<-droplevels(diff_data$CORR_RESP)
diff_data$STEP_PAIR_COLL<-droplevels(diff_data$STEP_PAIR_COLL)
diff_data$SUBJ_SHORT<-droplevels(diff_data$SUBJ_SHORT)
diff_data$TIME_LAG<-droplevels(diff_data$TIME_LAG)
diff_data$EXPERIMENT<-droplevels(diff_data$EXPERIMENT)


#make eight way forward difference coding
my.forward.diff = matrix(c(7/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,
                           6/8, 6/8,-2/8,-2/8,-2/8,-2/8,-2/8,-2/8,
                           5/8, 5/8, 5/8,-3/8,-3/8,-3/8,-3/8,-3/8,
                           4/8, 4/8, 4/8, 4/8,-4/8,-4/8,-4/8,-4/8,
                           3/8, 3/8, 3/8, 3/8, 3/8,-5/8,-5/8,-5/8,
                           2/8, 2/8, 2/8,2/8,2/8,2/8,-6/8,-6/8,
                           1/8,1/8,1/8,1/8,1/8,1/8,1/8,-7/8),ncol = 7)
contrasts(diff_data$STEP_PAIR_COLL)<-my.forward.diff

diff_data$TIME_LAG<-as.factor(diff_data$TIME_LAG)

my.forward.diff.4 = matrix(c(3/4, -1/4, -1/4, -1/4, 1/2, 1/2, -1/2, -1/2,
                             1/4, 1/4, 1/4, -3/4), ncol = 3)

contrasts(diff_data$TIME_LAG)<-my.forward.diff.4




#---------------------------------------

##### Frequentist model #####


library(lme4)

DIFF_step_pair_experiment<-glmer(diff_data$USER_CORR~diff_data$STEP_PAIR_COLL*diff_data$EXPERIMENT*diff_data$TIME_LAG
                                 +(1|diff_data$SUBJ_SHORT)+(1|diff_data$FILE1),
                                 family=binomial)

m <- DIFF_step_pair_experiment

print(summary(m),cor=F)
#does not converge, obviously.
#--------------------------------------

##### Bayesian model #####

# creating model matrices
x <- unname(model.matrix(~1+STEP_PAIR_COLL*EXPERIMENT+TIME_LAG, diff_data)) # matrix for fixed effects
attr(x, "assign") <- NULL
x_u <- unname(model.matrix(~1, data)) # matrix for random effects for subjects
attr(x_u, "assign") <- NULL
x_w <- unname(model.matrix(~1, data)) # matrix for items random effects
attr(x_w, "assign") <- NULL

# data list
stanDat <- list(accuracy = as.integer(diff_data$USER_CORR),         # dependent variable
                
                subj=as.numeric(factor(diff_data$SUBJ_SHORT)),  # subject id
                item=as.numeric(factor(diff_data$FILE1)),   # item id
                
                N_obs = nrow(data),                     # number of observations
                
                N_coef = ncol(x),                      # number of fixed effects
                N_coef_u = ncol(x_u),                    # number of random effects for subjects
                N_coef_w = ncol(x_w),                    # number of random effects for items
                
                x = x,                                 # fixed effects matrix
                x_u = x_u,                               # random effects matrix - subjects
                x_w = x_w,                               # random effects matrix - items
                
                N_subj=length(unique(diff_data$SUBJ_SHORT)),         # number of subjects
                N_item=length(unique(diff_data$FILE1)) )  # number of items

# model code
model_code_glmm <- "

data {
  int<lower=0> N_obs;                    //number of observations
  int<lower=0> N_coef;                   //fixed effects
  int<lower=0> N_coef_u;                 //random effects for subjects
  int<lower=0> N_coef_w;                 // random effects for items

  // subjects
  int<lower=1> subj[N_obs];          //subject id  
  int<lower=1> N_subj;               //number of subjects

  // items
  int<lower=1> item[N_obs];          //item id
  int<lower=1> N_item;               //number of items

  matrix[N_obs, N_coef] x;           //fixed effects design matrix
  matrix[N_obs, N_coef_u] x_u;         //subject random effects design matrix
  matrix[N_obs, N_coef_w] x_w;         //item random effects design matrix

  int accuracy[N_obs];              // accuracy
}

parameters {

  vector[N_coef] beta;                    // vector of fixed effects parameters 

  vector<lower=0> [N_coef_u] sigma_u;     // subject sd
  cholesky_factor_corr[N_coef_u] L_u;   // correlation matrix for random intercepts and slopes subj
  matrix[N_coef_u,N_subj] z_u;

  vector<lower=0> [N_coef_w] sigma_w;     // item sd
  cholesky_factor_corr[N_coef_w] L_w;    // correlation matrix for random intercepts and slopes item
  matrix[N_coef_w,N_item] z_w;

  real sigma_e;                // residual sd
}

transformed parameters{

  matrix[N_coef_u,N_subj] u;   // subjects random effects parameters 
  matrix[N_coef_w,N_item] w;   // items random effects parameters
  vector[N_obs] mu;

  // variance-covariance matrix for the random effects of subjects (intercept, slopes & correlations)
  {matrix[N_coef_u,N_coef_u] Lambda_u;
   Lambda_u = diag_pre_multiply(sigma_u, L_u);
   u = Lambda_u * z_u;
  }

  // variance-covariance matrix for the random effects of items (intercept, slopes & correlations)
  {matrix[N_coef_w,N_coef_w] Lambda_w;
   Lambda_w = diag_pre_multiply(sigma_w, L_w);
   w = Lambda_w * z_w;
  }

  mu = sigma_e + x * beta; // first define mu in terms of error (sigma_e) and fixed effects only (beta)

  for (i in 1:N_obs){
    for (uu in 1:N_coef_u)
      mu[i] = mu[i] + x_u[i,uu] * u[uu, subj[i]]; // adding to mu the subjects random effects part 
    for (ww in 1:N_coef_w)
      mu[i] = mu[i] + x_w[i,ww] * w[ww, item[i]]; // adding to mu the items random effects part
  }
}

model {

  // all priors here are weakly informative priors with a normal distribution (because the logit link function transforms the proportions into a normally distributed variable)
  // check what is the appropriate prior distribution for each of the variables in the data set and modify accordingly

  sigma_u ~ normal(0,1);
  sigma_w ~ normal(0,1);

  sigma_e ~ normal(0,10); // the mean (0) lies between -10 and 10 on logit scale, i.e., between very close to 0 and very close to 1 on probability scale (hence, it's a weakly informative prior)
  beta ~ normal(0,10);                    

  L_u ~ lkj_corr_cholesky(2.0); // 2.0 prior indicates no prior knowledge about the correlations in the random effects
  L_w ~ lkj_corr_cholesky(2.0);

  to_vector(z_u) ~ normal(0,1);
  to_vector(z_w) ~ normal(0,1);

  accuracy ~ bernoulli_logit(mu); // likelihood (the data)
}

// the following block generates some variables that might be interesting to look at in some cases, but not necessarily
generated quantities{

  matrix[N_coef_u,N_coef_u] Cor_u;
  matrix[N_coef_w,N_coef_w] Cor_w;

  int pred_correct[N_obs];
  real log_lik[N_obs];
  real diffP[N_coef];

  Cor_u = tcrossprod(L_u); // if you want to look at the subjects random effects correlations
  Cor_w = tcrossprod(L_w); // if you want to look at the items random effects correlations

  // the following loop translates the beta coefficients from logit scale (beta) to probability scale (diffP)
  for (j in 1:(N_coef)){
    diffP[j] = inv_logit(sigma_e + beta[j]) - inv_logit(sigma_e - beta[j]);
  }

  // generating the model's log likelihood to be used, for example, in model comparison
  for (i in 1:N_obs){
    pred_correct[i] = bernoulli_rng(inv_logit(mu[i]));
    log_lik[i] = bernoulli_logit_lpmf(accuracy[i]|mu[i]);
  }
}"

# fitting the model
library(rstan)
fit <- stan(model_code=model_code_glmm, 
            data=stanDat,
            iter=3000, # number of iterations in each chain
            chains=4, # number of chains
            control=list(adapt_delta=0.99, max_treedepth = 15) # this is not obligatory, only in order to facilitate model convergence and avoid divergent transitions
            )

# checking the chains
traceplot(fit, inc_warmup=F)

# checking that the Rhat is always 1 for all model parameters
model_sum <- summary(fit, probs=c(0.025,0.975))$summary
rhats <- model_sum[,"Rhat"]
round(summary(rhats),2) # everything should be 1 -- if not, model did not converge; in that case, try running the model again doubling the number of iterations per chain
( na <- model_sum[is.na(rhats),"Rhat"] ) # if parameters without Rhat are only those used in the random effects correlations, then it's not a problem

# saving the model output as a data frame (each column is an estimated posterior)
fit <- as.data.frame(fit)
# saving the data frame in R format
save(fit, file="fit.RData")




fit<-load(file="/Users/post-doc/Downloads/fit.RData")

fit<-as.data.frame(fit)

# taking only the fixed effects
beta <- fit[,grepl("beta",colnames(fit))]
      ###beta <- beta[,-1] # excluding the intercept


# creating a list of the fixed effects variables, in the order in which they were used in the model that was run
( cn <- colnames(model.matrix(~1+STEP_PAIR_COLL*EXPERIMENT+TIME_LAG, diff_data)) )
      ##cn <- cn[-1] # excluding the intercept

# creating an empty data frame that will be used in the plot
df <- data.frame(matrix(nrow=ncol(beta), ncol=0, data=NA, dimnames=list(c(),c())))
# creating a column "effect" with the parameter names
df$effect <- factor(cn)

# completing information in the data frame for plotting
for (i in 1:nrow(df)){
  
  # the mean of the posterior
  df[i,"mean"] <- mean(beta[,i])
  
  # probability that the posterior is smaller/greater than zero
  df[i,"probability_smaller"] <- mean(beta[,i]<0)
  df[i,"probability_bigger"] <- mean(beta[,i]>0)
  
  # range of the posterior (min / max)
  df[i,"min"] <- min(beta[,i])
  df[i,"max"] <- max(beta[,i])
  
  # 95% credible intervals
  df[i,"l95"] <- unname(quantile(beta[,i],probs=0.025))
  df[i,"h95"] <- unname(quantile(beta[,i],probs=0.975))
  
}

write.csv(df, file="beta_summary.csv")
library(ggplot2)
# plotting the posteriors
ggplot(data=df, aes(x=mean, y=effect)) + theme_bw() +
  geom_vline(aes(xintercept=0), size=1, linetype=2, col=gray(0.2)) + 
  geom_errorbarh(aes(xmax=max, xmin=min),height=0, size=1.3, col="#009E73") +
  geom_errorbarh(aes(xmax=l95, xmin=h95),linetype=1,height=0.2,size=1.5,col="#D55E00") +
  geom_point(size=2) + 
  scale_y_discrete(limits=rev(df$effect)) +
  theme(axis.title.y=element_text(size=12, angle=90),
        axis.title.x=element_text(size=24, angle=0),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=22, color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black")) + 
  ylab(" ") + xlab( expression(paste("Estimated difference (", hat(beta),")")))

# interpretation of the posteriors:
# since all variables were scaled around 0, 0 represents the point of "no effect".
# if the distribution does not include 0, then there's 100% probability that the effect is there.
# if the distribution includes 0, check whether it is included inside or outside the credible intervals
# check the probability that the posterior is greater/smaller than 0 to quantify the probability of the effect (i.e., to quantify the evidence for the effect)
