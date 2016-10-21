\name{stage_1_evaluation}
\alias{stage_1_evaluation}
\title{Stage 1-evaluation step of GSED}
\description{
\code{stage_1_evaluation} is used to evaluate the efficacy of the subgroup selected at the end of the first stage of GSED.
}
\usage{
stage_1_evaluation(keep, Z_1j, f, u)
}
\arguments{                
  \item{keep}{Vector containing the indices of the subgroups selected at stage 1.}
  \item{Z_1j}{Vector containing the Z-statistics (standard normal under H0) for each subgroup. Must be of length \code{N_subsets}.}
  \item{f}{Vector containing the prevalence rates of each subgroup.}
  \item{u}{Vector containing the upper boundaries for stagewise decisions.}
}
\value{
A list is returned, consisting of: 
\item{stage}{Integer containing the current step. Value is 1 by default, or -1 if the trial stops earlier for efficacy after this evaluation step.}
\item{S}{Vector containing the indices of the subgroups selected at stage 1 (=\code{keep}).}
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
stage_1_evaluation(keep=c(2,3), Z_1j=c(-0.49,1.07,1.44), f=c(0.6,0.2,0.2), u=c(2.7625,2.5204))
}
