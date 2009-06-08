\name{qls-internal}
\alias{cluster.size}
\alias{data.proc}
\alias{residual}
\alias{cormax.ar1}
\alias{cormax.exch}
\alias{cormax.markov}
\alias{cormax.tri}
\alias{cormax.fam}
\alias{cormax.exfam}
\alias{gen.zcor}
\alias{gee.ar1.fixed}
\alias{gee.exch.fixed}
\alias{gee.markov.fixed}
\alias{gee.tri.fixed}
\alias{gee.fam.fixed}
\alias{gee.exfam.fixed}
\alias{gee.fixed}
\alias{ar1.one}
\alias{ar1.two}
\alias{fam.one}
\alias{fam.two}
\alias{exch.one}
\alias{exch.two}
\alias{markov.one}
\alias{markov.two}
\alias{tri.one}
\alias{tri.two}
\alias{fam.one}
\alias{fam.two}
\alias{exfam.one}
\alias{exfam.two}
\alias{palp1}
\alias{palp2}
\alias{palp3}
\alias{palp4}
\alias{palp}
\title{Internal qls functions}
\description{
  Internal functions called by other functions in qlspack.
}
\usage{
cluster.size(id)
data.proc(data, formula, time = NULL, id, del.n)
residual(x, y, beta, family = "gaussian")
cormax.ar1(alpha, id, time = NULL)
cormax.exch(alpha, id, time = NULL)
cormax.markov(alpha, time)
cormax.tri(alpha, id, time = NULL)
cormax.fam(alpha, id, time = NULL)
cormax.exfam(alpha, id, time = NULL)
gen.zcor(cor.max, id, time = NULL, markov = FALSE)
gee.ar1.fixed(formula, data, id, alpha, family = "gaussian",
              time = NULL, std.err = std.err)
gee.exch.fixed(formula, data, id, alpha, family = "gaussian",
              time = NULL, std.err = std.err)
gee.markov.fixed(formula, data, id, alpha, family = "gaussian",
                 time = NULL, std.err = std.err)
gee.tri.fixed(formula, data, id, alpha, family = "gaussian",
              time = NULL, std.err = std.err)
gee.tri.fixed(formula, data, id, alpha, family = "gaussian",
              time = NULL, std.err = std.err)
gee.exfam.fixed(formula, data, id, alpha, family = "gaussian",
                time = NULL, std.err = std.err)
gee.fixed(formula, data, id, alpha, family = "gaussian", time = NULL,
          correlation, std.err = std.err)
ar1.one(resid, time = NULL, id)
ar1.two(alpha, resid=NULL, time = NULL, id = NULL)
exch.one(resid, time = NULL, id)
exch.two(alpha, resid, time = NULL, id)
markov.one(resid, time = NULL, id)
markov.two(alpha, resid, time = NULL, id)
tri.one(resid, time = NULL, id)
tri.two(alpha, resid, time = NULL, id)
exfam.one(resid, time = NULL, id)
exfam.two(alpha, id, time = NULL)
exfam.one(resid, time = NULL, id)
exfam.two(alpha, id, time = NULL)
palp1(t)
palp2(t)
palp3(t)
palp4(t)
palp(t, j) }
\details{
 These are not to be called directly by the user.
}
\keyword{internal}
