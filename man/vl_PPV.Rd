% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vl_PPV.R
\name{vl_PPV}
\alias{vl_PPV}
\title{Positive Predicted Value curve}
\usage{
vl_PPV(
  predicted,
  label,
  Nleft = 100,
  plot = FALSE,
  xlim = NULL,
  ylim = NULL,
  xlab = "Prediction score",
  ylab = "Positive pred. value (\%)",
  lty.1 = 1,
  lty.2 = 3,
  col = "black",
  cex.max = 0.7,
  pos.max = 3,
  offset.max = 0.5,
  add = FALSE,
  ...
)
}
\arguments{
\item{predicted}{Predicted values from the model (ranging from 0 to 1).}

\item{label}{A vector of logical labels (or that can be coerced to logical).}

\item{Nleft}{Number of enhancers left before cutoff. Default= 100.}

\item{plot}{Should the PPV be plotted?}

\item{xlim}{x limits for plotting. Default= NULL.}

\item{ylim}{y limits for plotting. Default= NULL.}

\item{xlab}{x label. Default= "Prediction score".}

\item{ylab}{y label. Default= "Positive pred. value (percentage)".}

\item{lty.1}{Line type before Nleft. Default= 1.}

\item{lty.2}{Line type after Nleft. Default= 3.}

\item{col}{Color.}

\item{cex.max}{cex for the max PPV text label.}

\item{pos.max}{pos for the PPV text label Default= 3.}

\item{offset.max}{Offset for plotting the PPV text label. Default= .5.}

\item{add}{When plot is set to TRUE, should only the lines/label be added to an existing plot?}

\item{...}{Extra arguments to be passed to lines}
}
\value{
PPV plot
}
\description{
Positive Predicted Value curve
}
