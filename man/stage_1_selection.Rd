\name{stage_1_selection}
\alias{stage_1_selection}
\title{Stage 1-selection step of GSED}
\description{
\code{stage_1_selection} is used to determine the subgroup selected at the end of the first stage of GSED.
}
\usage{
stage_1_selection(N_subsets, Z_1j, l, ordering, increasing_theta=FALSE)
}
\arguments{                
  \item{N_subsets}{Integer representing the number of possible subgroups.}
  \item{Z_1j}{Vector containing the Z-statistics (standard normal under H0) for each subgroup. Must be of length \code{N_subsets}.}
  \item{l}{Vector containing the lower boundaries for stagewise decisions.}
  \item{ordering}{Boolean indicating if the subgroups (theta) are ordered.}
  \item{increasing_theta}{Boolean indicating if greater values of theta parameters represent better treatment effects. The default value is set at FALSE.}
}
\value{
A vector containing the indices of the subgroups selected is returned.
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
stage_1_selection(N_subsets=3, Z_1j=c(-0.49,1.07,1.44), l=c(0.7962,2.5204), ordering=FALSE)
}
