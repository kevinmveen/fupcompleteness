\name{fupcompleteness}
\alias{hello}
\title{Follow-up completness}
\usage{
hello()
}
\description{
Prints 'Returns a list with follow-up completeness and data frame to calculate it'.
}
\Arguments{test}
\examples{
dat.output2 = fupcompletness::fup.completeness(date.inclusion = input_data$date.inclusion,
                                 last.fup.date = input_data$last.fup.date,
                                 end.date = input_data$end.date,
                                 status = input_data$status,
                                 cencode = 0, 
                                 method = "mclarkc")
}
