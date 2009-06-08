\name{rat}
\alias{rat}
\docType{data}
\title{ ~~ data name/kind ... ~~}
\description{
  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(rat)}
\format{
  A data frame with 291 observations on the following 6 variables.
  \describe{
    \item{\code{id2}}{a numeric vector}
    \item{\code{id}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{group}}{a numeric vector}
    \item{\code{bp}}{a numeric vector}
    \item{\code{highbp}}{a numeric vector}
  }
}
\details{
The description of each column of the data set is as follows:

id2 = the id variable for each rat that is provided in Davis (2002).

id = a new id variable that takes value 1,2,..43 after sorting on id and group.

time = the timing of each measurement.

group = the group variable that takes value 1, 2, 3, or 4

bp = the blood pressure value.

highbp = a variable that takes value 1 if the rat's blood pressure is at least 100.
}
\source{
Table 6.11 of Davis (2002)
}
\references{
Davis, C.(2002). \emph{Statistical Methods for the Analysis of Repeated Measurements.}
}
\examples{
data(rat)
## maybe str(rat) ; plot(rat) ...
}
\keyword{datasets}
