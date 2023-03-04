\name{fup.completeness}
\alias{fup.completeness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to calculate follow-up completness in (medical) cohort studies
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to caculate follow-up complentess according to the Clark C (1), Modified Clark C (2), percentage method (3), follow-up index (4), Person-Time Follow-up Rate (FPT) (3) and simplified person-time method (SPT) (3)
}
\usage{
fup.completeness(date.inclusion, end.date, last.fup.date, status, cencode =0, death.date = NULL,            death = NULL, deathcode = 1, method = "clarkc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date.inclusion}{Vector of dates of inclusion of the subjects. Should be as.Date("2000-01-01") format.}
  \item{end.date}{Vector of end dates of the study. Should be as.Date("2000-01-01") format.}
  \item{last.fup.date}{Vector of date of last follow-up of the subjects. Cannot be later than end date of the study. Should be as.Date("2000-01-01") format.}
  \item{status}{Vector of status of subjects at last follow-up date. Usually 1 = event of interest and 0 = censoring}
  \item{cencode}{Code used for censoring in the status vector. Default = 0}
  \item{death.date}{Vector of death dates of the subjects. If not entered status is assumed to be death. Default =NULL}
  \item{deathcode}{Code for death in status vector. Default = 1}
  \item{method}{Method used to caculate follow-up completness. This can be "clarkc" (Clark C), "mclarkc" (modified Clark C), "percent" (percentage method), "FUI" (follow-up index), "PTFR", (Person-Time Follow-up Rate ), "SPT" (simplified person-time method)}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List containing:
\item{dat}{Data frame of inputted data, with extra colums needed for the caculations}
\item{fupc}{Numerical of follow-up completness}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin M Veen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(simsurv)
library(dplyr)
library(survival)

# Define data
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

# add random censoring
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

last.fup.date = rep(as.Date("2007-01-01"),1000) + dat2$time*365.25

# define dates
input_data = data.frame(date.inclusion = rep(as.Date("2007-01-01"),1000),
                        end.date = rep(as.Date("2012-01-01"),1000),
                        last.fup.date = rep(as.Date("2007-01-01"),1000) + dat2$time*365.25,
                        status = dat2$status2)
#run function
dat.output2 = fupcompletness::fup.completeness(date.inclusion = input_data$date.inclusion,
                                 last.fup.date = input_data$last.fup.date,
                                 end.date = input_data$end.date,
                                 status = input_data$status,
                                 cencode = 0,
                                 method = "mclarkc")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
