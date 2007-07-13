\name{qls-internal}
\alias{cluster.size}
\alias{proc}
\alias{data.process}
\alias{data.proc.exfam}
\alias{residual}
\alias{cormax.ar1}
\alias{cormax.exch}
\alias{cormax.markov}
\alias{cormax.tri}
\alias{cormax.exfam}
\alias{gen.zcor}
\alias{gee.ar1.fixed}
\alias{gee.exch.fixed}
\alias{gee.markov.fixed}
\alias{gee.tri.fixed}
\alias{gee.exfam.fixed}
\alias{gee.fixed}
\alias{ar1.one}
\alias{ar1.two}
\alias{exch.one}
\alias{exch.two}
\alias{markov.one}
\alias{markov.two}
\alias{tri.one}
\alias{tri.two}
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
proc(data = "NA", time = "NA", id)
data.process(data = "NA", time = "NA", id)
data.proc.exfam(data = "NA", time = "NA", id)
residual(x, y, beta, family = "gaussian")
cormax.ar1(alpha, id, time = "NA")
cormax.exch(alpha, id, time = "NA")
cormax.markov(alpha, time)
cormax.tri(alpha, id, time = "NA")
cormax.exfam(alpha, id, time = "NA")
gen.zcor(cor.max, id, time = "NA", markov = FALSE)
gee.ar1.fixed(formula, data, id, alpha, family = "gaussian",
              time = "NA", std.err = std.err)
gee.exch.fixed(formula, data, id, alpha, family = "gaussian",
              time = "NA", std.err = std.err)
gee.markov.fixed(formula, data, id, alpha, family = "gaussian",
                 time = "NA", std.err = std.err)
gee.tri.fixed(formula, data, id, alpha, family = "gaussian",
              time = "NA", std.err = std.err)
gee.exfam.fixed(formula, data, id, alpha, family = "gaussian",
                time = "NA", std.err = std.err)
gee.fixed(formula, data, id, alpha, family = "gaussian", time = "NA",
          correlation, std.err = std.err)
ar1.one(resid, time = "NA", id)
ar1.two(alpha, resid = "NA", time = "NA", id = "NA")
exch.one(resid, time = "NA", id)
exch.two(alpha, resid = "NA", time = "NA", id)
markov.one(resid, time = "NA", id)
markov.two(alpha, resid = "NA", time = "NA", id)
tri.one(resid, time = "NA", id)
tri.two(alpha, resid = "NA", time = "NA", id)
exfam.one(resid, time = "NA", id)
exfam.two(alpha, id, time = "NA")
palp1(t)
palp2(t)
palp3(t)
palp4(t)
palp(t, j) }
\details{
 These are not to be called directly by the user.
}
\keyword{internal}
