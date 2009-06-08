\name{qls}
\alias{qls}
\title{Fit Quasi Least Squares (QLS) Estimating Equations}
\description{
  The qls function fits quasi least square estimating
  equations based on the geeglm function in the geepack
  and cor.estimate funcion in the qlspack. qls has a syntax
  similar to glm and returns an object similar to a glm
  object. An important feature of qls, is that an anova method
  exists for these models.
}
\usage{
qls(formula, data, id, family = "gaussian",
time = NULL, correlation = "ar1", std.err = "san.se")
}
\arguments{
  \item{formula}{
The model to be fitted. The form is similar to the item
documentated in \code{geeglm}. 
          }
  \item{data}{
A data frame containing the variables in the model.
             }
  \item{id}{
a vector which identifies the clusters.  The length of `id'
should be the same as the number of observations.  Data are
assumed to be sorted so that observations on a cluster are
contiguous rows for all entities in the formula. The 'id's
for different clusters should be different, but need not to
be consecutive.
           }
  \item{family}{
A character string describing the error distribution and link
function to be used in the model. There are three options:
"guassian", "binomial" and "poisson". The default option
is "guassian".
               }
  \item{time}{
a vector which identifies the time in the clusters. The length
of 'time' should be the same as the number of observations.
This argument is used if and only if 'correlation == "markov"'.
             }

  \item{correlation}{
a character string specifying the correlation structure. The
following are permitted: '"ar1"', '"exchangeable"', '"markov"',
'"tridiagonal"', '"fam"' and '"ex.fam"'.
                     }

  \item{std.err}{
See corresponding documentation to \code{geeglm}.
                 }}
\value{
  An object of type 'qlsglm'.
}

\references{
Chaganty, N. R. 1997. An alternative approach to the
analysis of longitudinal data via generalized estimating equations.
\emph{Journal of Statistical Planning and Inference} \bold{63}: 39--54.

Shults, J. 1996. The analysis of unbalanced and unequally spaced
longitudinal data using quasi-least squares.
Ph.D. Thesis, Department of Mathematics and Statistics,
Old Dominion University: Norfolk, Virginia.

Shults, J. and Chaganty, N.R. 1998. Analysis of serially
correlated data using quasi-least squares.
\emph{Biometrics} \bold{54}: 1622--1630.

Chaganty, N.R. and Shults, J. 1999. On eliminating the asymptotic
bias in the quasi-least squares estimate of the correlation parameter.
\emph{Journal of Statistical Planning and Inference} \bold{76}: 127--144.
}
\author{Jichun Xie, jichun@mail.med.upenn.edu}
\note{qls only works for complete data. Thus if there are NA's in data you can specify data=na.omit(mydata).
}

\section{Warning }{qls has not been thoroughly tested. Please report bugs.}

\seealso{\code{\link{glm}}}
\examples{
require(qlspack)
data(rat)
qlsfit.fam <- qls(bp ~ time + as.factor(group), data = rat, id = rat$id,
          time = rat$time, correlation = "fam")
summary(qlsfit.fam)
}
\keyword{models}
