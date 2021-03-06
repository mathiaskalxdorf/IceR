% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General.R
\name{load_MaxQ_data}
\alias{load_MaxQ_data}
\title{Load MaxQuant result files}
\usage{
load_MaxQ_data(
  path = NA,
  min_pep_count = 1,
  min_pep_count_criteria = c("all", "unique"),
  remove_contaminants = T,
  remove_reverse = T,
  intensity_used = c("LFQ intensity", "iBAQ", "Intensity")
)
}
\arguments{
\item{path}{Optional path to folder containing MaxQuant outputs. By default set to NA. In this case a file browser is opened. If path is directly specified, it has to end with \\}

\item{min_pep_count}{Minimal required number of quantified peptides per protein. By default set to 1.}

\item{min_pep_count_criteria}{Criteria how to count quantified peptides per protein. Either all or only unique peptides are counted. By default set to "all".}

\item{remove_contaminants}{Boolean value indicating if contaminants should be removed. By default set to T.}

\item{remove_reverse}{Boolean value indicating if reverse hits should be removed. By default set to T.}

\item{intensity_used}{Specifying which protein quantification data should be used. Selection between "LFQ intensity", "iBAQ" or "Intensity". By default set to "LFQ intensity". Requires corresponding quantification results to be calculated by MaxQuant and stored in respective columns.}
}
\value{
List object containing protein and peptide quantification information in sub-lists named Protein_level and Peptide_level, respectively.
}
\description{
Load MaxQuant result files
}
\details{
Wrapper function to load and filter MaxQuant results.
}
