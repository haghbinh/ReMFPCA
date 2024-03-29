% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvmfdClass.R
\name{mvmfd}
\alias{mvmfd}
\alias{Mvmfd}
\title{Define a Set of Multivariate Multidimensional Functional Data objects}
\usage{
Mvmfd(...)
}
\arguments{
\item{...}{A `mfd` objects which have separated by comma}
}
\description{
The `mvmfd` class represents functional data ...

Constructor for `mvmfd` objects (same as `Mvmfd`)
}
\examples{
require(fda)
bs1 <- create.fourier.basis(c(0, 2 * pi), 5)
bs2 <- create.bspline.basis(c(0, 1), 7)
bs3 <- create.exponential.basis(c(0, 2), 3)
nobs <- 10
argval1 <- seq(0, 2 * pi, length.out = 12)
X1 <- outer(sin(argval1), seq(0.5, 1.5, length.out = nobs))
mdbs1 <- Basismfd(bs1)
mfd1 <- Mfd(argval1, X1, mdbs1)
mdbs2 <- Basismfd(bs1)
argval2 <- argval1
X2 <- outer(cos(argval2), seq(0.2, 1.5, length.out = nobs))
mfd2 <- Mfd(argval2, X2, mdbs1)
mvmfd1 <- Mvmfd(mfd1, mfd2)
mvmfd1[1]
mvmfd1[1, 1]
mvmfd1[1:5, 2]
mvmfd1[, 1]
mvmfd1[1:5, ]
evalarg <- list(argval1, argval2)
mvmfd1$eval(evalarg)
mvmfd1 + mvmfd1
mean(mvmfd1)
inprod_mvmfd(mvmfd1, mvmfd1)
norm_mvmfd(mvmfd1)
plot(mvmfd1)
bimfdplot(mvmfd1)

}
\seealso{
\code{\link{mvbasismfd}}, \code{\link{mfd}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{basis}}{A `mvbasismfd` object}

\item{\code{coefs}}{a matrix of the coefficients.}

\item{\code{nobs}}{number of observation}

\item{\code{nvar}}{number of variables}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-mvmfd-new}{\code{mvmfd$new()}}
\item \href{#method-mvmfd-eval}{\code{mvmfd$eval()}}
\item \href{#method-mvmfd-print}{\code{mvmfd$print()}}
\item \href{#method-mvmfd-clone}{\code{mvmfd$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mvmfd-new"></a>}}
\if{latex}{\out{\hypertarget{method-mvmfd-new}{}}}
\subsection{Method \code{new()}}{
Constructor for `mvmfd` objects (same as 'Mvmfd')
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mvmfd$new(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{A `mfd` objects which have separated by comma}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mvmfd-eval"></a>}}
\if{latex}{\out{\hypertarget{method-mvmfd-eval}{}}}
\subsection{Method \code{eval()}}{
Eval method for `mvmfd` objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mvmfd$eval(evalarg)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{evalarg}}{A list of numeric vectors of argument values at which the `mvmfd` is to be evaluated.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of evaluated values
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mvmfd-print"></a>}}
\if{latex}{\out{\hypertarget{method-mvmfd-print}{}}}
\subsection{Method \code{print()}}{
Print method for `mvmfd` objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mvmfd$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Additional arguments to be passed to `print`}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mvmfd-clone"></a>}}
\if{latex}{\out{\hypertarget{method-mvmfd-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mvmfd$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
