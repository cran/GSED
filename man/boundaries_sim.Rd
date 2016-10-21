\name{boundaries_sim}
\alias{boundaries_sim}
\title{Lower and upper boundaries for GSED}
\description{
\code{boundaries_sim} is used to estimate lower and upper boundaries for GSED based on simulations of trials.
}
\usage{
boundaries_sim(K_stages, N_subsets, f, ratio_Delta_star_d1, ordering, 
increasing_theta=FALSE, seed=42, n_trials, alpha_spending, 
one_minus_alpha_spending)
}
\arguments{                
  \item{K_stages}{Integer indicating the number of stages in the design.}
  \item{N_subsets}{Integer representing the number of possible subgroups.}
  \item{f}{Vector containing the prevalence rates of each subgroup. Must be of length \code{N_subsets}.}
  \item{ratio_Delta_star_d1}{Vector containing the ratio between the (observed Fisher) information increments at each stage >1 with the (observed Fisher) information at stage 1. Must be of length \code{K_stages}-1.}
  \item{ordering}{Boolean indicating if the subgroups (theta) are ordered.}
  \item{increasing_theta}{Boolean indicating if greater values of theta parameters represent better treatment effects. The default value is set at FALSE.}
  \item{seed}{Interger representing the seed. The default value is set at 42.}
  \item{n_trials}{Integer indicating the number of trials to simulate.}
  \item{alpha_spending}{Vector containing the values of the alpha-spending function at each time of the analysis (including 0 at time 0 and alpha at time 1). Must be of length \code{K_stages}+1.}
  \item{one_minus_alpha_spending}{Vector containing the values of the 1-alpha-spending function at each time of the analysis (including 0 at time 0 and 1-alpha at time 1). Must be of length \code{K_stages}+1.}
}
\value{
A list is returned, consisting of two vectors containing the lower and upper boundaries:
\item{l}{Vector of lower boundaries at each stage.}
\item{u}{Vector of upper boundaries at each stage.}
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
#For testing purpose only, larger number of simulations required (see in comments below)
boundaries_sim(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), ratio_Delta_star_d1=c(1), 
ordering=FALSE, seed=42, n_trials=3, alpha_spending=c(0,0.0125,0.025), 
one_minus_alpha_spending=c(0,0.4875,0.975))

#boundaries_sim(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), ratio_Delta_star_d1=c(1), 
#ordering=FALSE, seed=42, n_trials=10000000, alpha_spending=c(0,0.0125,0.025), 
#one_minus_alpha_spending=c(0,0.4875,0.975))
}
