#simulate data
rm(list = ls())
library(simsurv)
library(dplyr)
library(survival)
library(survminer)



cov <- data.frame(id = 1:1000,
                  trt = rbinom(1000, 1, 0.5))

# Simulate the event times
dat <- simsurv(lambdas = 0.1,
               gammas = 1.5,
               betas = c(trt = -0.5),
               x = cov,
               maxt = 5)

# Merge the simulated event times onto covariate data frame
dat <- merge(cov, dat)


head(dat)

KM = survfit(Surv(time = dat$eventtime, event = dat$status ==1) ~ 1, data = dat)
ggsurvplot(KM, data = dat, conf.int = T)

dat2 = dat
dat2$cens = rbinom(1000, 1, 0.3)
dat2$time = dat2$eventtime

for(i in 1:nrow(dat2)){
  if(dat2$cens[i] == 1) {
    dat2$time[i] = runif(1, min = 1/365.25, max = dat2$eventtime[i])
  }
}

dat2$status2 = dat2$status
dat2$status2 = replace(dat2$status2, dat2$status == 1 & dat2$eventtime > dat2$time, 0)

dat$lossfup = "No loss to FUP"
dat2$lossfup = "Loss to FUP (30% of patients)"

dat3 = cbind.data.frame(event = c(dat$status, dat2$status2),
                        time = c(dat$eventtime, dat2$time),
                        lossfup = c(dat$lossfup, dat2$lossfup))

head(dat3)
KM = survfit(Surv(time = dat3$time, event = dat3$event ==1) ~ lossfup, data = dat3)
ggsurvplot(KM, data = dat3, conf.int = T, pval = T)


coxph(Surv(time = dat$eventtime, event = dat$status ==1) ~ trt, data = dat)
coxph(Surv(time = dat2$time, event = dat2$status2 ==1) ~ trt, data = dat2)
2007+5

last.fup.date = rep(as.Date("2007-01-01"),1000) + dat2$time*365.25

input_data = data.frame(date.inclusion = rep(as.Date("2007-01-01"),1000),
                        end.date = rep(as.Date("2012-01-01"),1000),
                        last.fup.date = rep(as.Date("2007-01-01"),1000) + dat2$time*365.25,
                        status = dat2$status2)




