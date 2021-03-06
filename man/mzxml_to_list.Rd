% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature_alignment_requantification.R
\name{mzxml_to_list}
\alias{mzxml_to_list}
\title{Prepare mzXML files for IceR workflow}
\usage{
mzxml_to_list(path_to_mzXML, n_cores = 2)
}
\arguments{
\item{path_to_mzXML}{Path to folder containing mzXML files}

\item{n_cores}{Numbers of CPU cores which should be used to perform conversion}
}
\value{
Resulting ion tables are stored in a sub-directory (all_ion_lists) of the mzXML folder as .RData files
}
\description{
Prepare mzXML files for IceR workflow
}
\details{
Converts MS1 spectra in mzXML files into tables containing m/z, RT and intensity information per ion
}
