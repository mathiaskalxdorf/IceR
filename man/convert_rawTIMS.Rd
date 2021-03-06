% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature_alignment_requantification.R
\name{convert_rawTIMS}
\alias{convert_rawTIMS}
\title{Extract MS1-spectra from raw TIMS-ToF Pro data}
\usage{
convert_rawTIMS(path_to_raw = NULL)
}
\arguments{
\item{path_to_raw}{Path to folder containing d-files}
}
\value{
Resulting ion tables are stored in a sub-directory (all_ion_lists) of the raw folder as .RData files
}
\description{
Extract MS1-spectra from raw TIMS-ToF Pro data
}
\details{
Convert MS1 spectra with TIMS information in raw files of Bruker TIMS-ToF Pro Mass Spectrometers into tables containing m/z, RT, inverse ion mobility and intensity information per ion. This process is currently very slow and can take several days even for a small data set. It also requires enough space on the home drive of the PC as well as at least 50 Gb of memory.
}
