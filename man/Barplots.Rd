% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General.R
\name{Barplots}
\alias{Barplots}
\title{Generate barplots}
\usage{
Barplots(
  Data,
  ErrbarData = NA,
  Name = "",
  xlab = "X-axis",
  ylab = "Y-axis",
  main = "Titel",
  col = "lightblue",
  AvgLine = T,
  digits_average = 0,
  Legends = NA,
  Legendscol = NA,
  Legendtitle = "Legend",
  Legendpos = "topright",
  shownumbers = T,
  shownumbers_digits = 1,
  ylim = NA,
  logy = F,
  margins = c(10.1, 4.1, 4.1, 4.1),
  inset = c(-0.1, 0)
)
}
\arguments{
\item{Data}{Numeric vector or table of samples in columns}

\item{ErrbarData}{Data for errorbars}

\item{Name}{Names}

\item{xlab}{X-Axis label}

\item{ylab}{Y-Axis label}

\item{main}{Plot main title}

\item{col}{Color}

\item{AvgLine}{Show average line?}

\item{digits_average}{Number of digits of average indication}

\item{Legends}{Legends}

\item{Legendscol}{Color of legends}

\item{Legendtitle}{Titel for legends}

\item{Legendpos}{Legend position}

\item{shownumbers}{Show numbers on top of bars}

\item{shownumbers_digits}{Number of digits for shown numbers}

\item{ylim}{y-axis limits}

\item{logy}{Y-Axis in log-scale?}

\item{margins}{Margins}

\item{inset}{Inset for legend}
}
\value{
Plot.
}
\description{
Generate barplots
}
\details{
Generate barplots
}
