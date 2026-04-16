#Inferencia 2026
#R Codes


require(truncnorm)
require(tidyverse)
require(car)
dat=read.csv("C:\\Users\\SOURAV\\Downloads\\Hackathon 2026\\trialx_hackathon_data.csv")
dat$sex=as.factor(dat$sex)
dat$site=as.factor(dat$site)
dat$treatment=as.factor(dat$treatment)
dat$drop_flag=ifelse(is.na(dat$dropout_week),0,1)

#Scrutiny of Data
dat[which(is.na(dat$age)),]

dat[which(is.na(dat$bmi)),]


#Distributional Checking of Key Variables
par(mfrow=c(1,2))

hist(dat$bmi,prob=T,main = "Histogram of BMI",
xlab="BMI measure",
     ylab="Proportion of Patients")
curve(dnorm(x,mean=mean(dat$bmi,na.rm=T),sd=sd(dat$bmi,na.rm=T)),add=T)

hist(dat$age,prob=T,main = "Histogram of Age",
xlab="Age of Patients",
     ylab="Proportion of Patients")
curve(dnorm(x,mean=mean(dat$age,na.rm=T),sd=sd(dat$age,na.rm=T)),add=T)


#Distributional Checking of Key Variables (Formal Test)
shapiro.test(dat$bmi)

shapiro.test(dat$age)


#Variable Pre-processing
impute=function(x){
  mu=mean(x,na.rm=TRUE)
  std=sd(x,na.rm=TRUE)
  low=min(x,na.rm=TRUE)
  up=max(x,na.rm=TRUE)
  #identify NA indices
  na_id=is.na(x)
  #Generate values from Truncated Normal Distribution
  x[na_id]=rtruncnorm(sum(na_id),a=low,b=up,mean=mu,sd=std)
  return(x)
}


dat$age=round(impute(dat$age),0)

dat$bmi=round(impute(dat$bmi),2)


#EDA
#Data Quality Report
tukey_outliers <- function(df, column) {
  values <- na.omit(df[[column]])
  
  q1 <- quantile(values, 0.25, na.rm = TRUE)
  q3 <- quantile(values, 0.75, na.rm = TRUE)
  iqr <- (q3 - q1)
  
  lower_bound <- q1 - (1.5 * iqr)
  upper_bound <- q3 + (1.5 * iqr)
  
  outliers <- values[values < lower_bound | values > upper_bound]
  
  return(outliers)
}


v <- list()

for(i in names(dat)){
  if(is.numeric(dat[[i]])){
    out <- tukey_outliers(dat, i)
    
    if(length(out) > 0){
      v[[i]] <- out
      cat("\nOutliers in", i, ":\n")
      print(out)
    }
  }
}


#outliers via boxplot
vars <- c("age", "bmi", "sbp_0", "sbp_4", "sbp_8", "sbp_12")

boxplot(dat[ , vars],
main = "Boxplots with Outliers", cex.main=0.5,
xlab = "Variables",
ylab = "Values", cex=0.5,
col = "lightblue", cex.lab=0.5, cex.axis=0.5,
border = "black",
las = 2)



#checking the Balance of baseline
#covariates (age, sex, BMI) across treatment arms.
require(lmtest)
require(naniar)
require(tableone)
require(sandwich)
require(multcomp)

#Table to check balance of baseline covariates
vars <- c("age", "sex", "bmi")
table1 <- CreateTableOne(vars = vars, strata = "treatment", data = dat)
print(table1, showAllLevels = TRUE)


#protocol deviation recorded at site 3
table(dat$site,dat$deviation_flag)


#Trt across gender
tab=table(dat$sex,dat$treatment)

barplot(tab,beside=TRUE,main="Assignment of Treatment across Gender",xlab = "Gender (Violet-Female,Yellow-Male)",
        ylab="Number of patients assigned",col=c("violet","yellow"))



#sex and treatment assignment were independent?
tab=table(dat$sex,dat$treatment)
tab
chisq.test(tab)



#Distributional checks per site and treatment arm
library(ggplot2)
library(patchwork)


g1=ggplot(dat, aes(x = age, fill = factor(site))) +
geom_density(alpha = 0.4)

g2=ggplot(dat, aes(x = bmi, fill = factor(site))) +
geom_density(alpha = 0.4)

g3=ggplot(dat, aes(x = age, fill = factor(treatment))) +
geom_density(alpha = 0.4)

g4=ggplot(dat, aes(x = bmi, fill = factor(treatment))) +
geom_density(alpha = 0.4)

(g1 | g2) / (g3 | g4)



#Inferential Analysis
#Treatment Effect Definition

dat$del_sbp=(dat$sbp_0-dat$sbp_12)

#Normality Check
dat$del_sbp=(dat$sbp_0-dat$sbp_12)
shapiro.test(dat$del_sbp)


#Model Specification
fit=lm(del_sbp~sex+site+age+bmi+treatment+sbp_0,data=dat)

coeftest(fit, vcov = vcovHC(fit))



#Device Switch at Site - 2
a=(dat[which(dat$site==2),]$sbp_0+dat[which(dat$site==2),]$sbp_4)/2;a=na.omit(a)
b=(dat[which(dat$site==2),]$sbp_8+dat[which(dat$site==2),]$sbp_12)/2;b=na.omit(b)

#Normallity check
shapiro.test(a);shapiro.test(b)


var.test(a,b)

t.test(a,b,var.equal = TRUE)


#Dropout Mechanism
drop_model <- glm(drop_flag ~ age + sex + bmi + sbp_0 + treatment + site,
data = dat,
family = binomial)
summary(drop_model)

round(exp(coef(drop_model)),3)


#Site-level Analysis with Multiplicity Correction
site_results <- dat %>%
group_by(site) %>%
do(model = lm(del_sbp ~ treatment + sbp_0 + age + sex + bmi, data = .))
pvals <- sapply(site_results$model, function(m) summary(m)$coeff["treatmentTrX", "Pr(>|t|)"])
p.adjust(pvals, method = "bonferroni")


#Primary Hypothesis Test
robust <- coeftest(fit, vcov = vcovHC(fit))
#Extract TrX row
trx <- robust["treatmentTrX", ]
#Degrees of freedom
df_res <- df.residual(fit)
#CI
ci <- coefci(fit, vcov = vcovHC(fit))["treatmentTrX", ]
#Effect size (Cohen's d approximation)
#d = beta / residual SD
effect_size <- coef(fit)["treatmentTrX"] / sigma(fit)
#Create table
result_table <- data.frame(
Estimate = trx["Estimate"],
Robust_SE = trx["Std. Error"],
t_value = trx["t value"],
df = df_res,
p_value = trx["Pr(>|t|)"],
CI_lower = ci[1],
CI_upper = ci[2],
Effect_Size_d = effect_size
)
round(result_table,6)


#Robustness and Sensitivity
#Sensitivity Analysis Excluding Site 3

df_sens <- dat %>% filter(site != 3)
fit_sens <- lm(del_sbp ~ treatment + sbp_0 + age + sex + bmi + factor(site),
data = df_sens)
summary(fit_sens)

orig <- coef(fit)["treatmentTrX"]
sens <- coef(fit_sens)["treatmentTrX"]
abs_change <- abs(sens - orig)
perc_change <- (abs_change / abs(orig)) * 100

abs_change
perc_change















































































































































































































































































