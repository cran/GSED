\name{max_FI}
\alias{max_FI}
\title{Maximum Fisher information}
\description{
\code{max_FI} is used to estimate maximum Fisher information based on two power criteria. 
- The first criterion consider the maxmimum Fisher information such that there is a pre-defined power to declare efficacy in the entire population for a given vector of parameters representing treatment effetcs in each subgroup.
- The second criterion consider the maxmimum Fisher information such that there is a pre-defined power to declare efficacy in at least one subgroup for a given vector of parameters representing treatment effetcs in each subgroup.
}
\usage{
max_FI(K_stages, N_subsets, f, ratio_Delta_star_d1, l, u, theta, pow, ordering, 
increasing_theta=FALSE, seed=42, n_trials, rule)
}
\arguments{ 
  \item{K_stages}{Integer indicating the number of stages in the design.}
  \item{N_subsets}{Integer representing the number of possible subgroups.}
  \item{f}{Vector containing the prevalence rates of each subgroup. Must be of length \code{N_subsets}.}
  \item{ratio_Delta_star_d1}{Vector containing the ratio between the (observed Fisher) information increments at each stage >1 with the (observed Fisher) information at stage 1. Must be of length \code{K_stages}-1.}
  \item{l}{Vector containing the lower boundaries for stagewise decisions. Must be of length \code{K_stages}.}
  \item{u}{Vector containing the upper boundaries for stagewise decisions. Must be of length \code{K_stages}.}
  \item{theta}{Vector of parameters representing treatment effects in each subgroup. Must satisfy the properties detailed in Magnusson and Turnbull's article (reparametrization can be needed).}
  \item{pow}{Value representing the desired power.}
  \item{ordering}{Boolean indicating if the subgroups (theta) are ordered.}
  \item{increasing_theta}{Boolean indicating if greater values of theta parameters represent better treatment effects. The default value is set at FALSE.}
  \item{seed}{Interger representing the seed. The default value is set at 42.}
  \item{n_trials}{Integer indicating the number of trials to simulate.}
  \item{rule}{Integer with value either 1 or 2 for power criteria detailed in description section (1 for entire population, 2 for at least one subgroup).}
}
\value{
A value representing the maximum Fisher information is returned.
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
#For testing purpose only, larger number of simulations required (see in comments below)
max_FI(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), ratio_Delta_star_d1=c(1), l=c(0.7962, 2.5204), 
u=c(2.7625, 2.5204), theta=c(0.2,0.2,0.2), pow=0.9, ordering=FALSE, increasing_theta=FALSE, 
seed=140691, n_trials=3, rule=1)
       
#max_FI(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), ratio_Delta_star_d1=c(1), l=c(0.7962, 2.5204), 
#u=c(2.7625, 2.5204), theta=c(0.2,0.2,0.2), pow=0.9, ordering=FALSE, increasing_theta=FALSE, 
#seed=140691, n_trials=10000000, rule=1)

#max_FI(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), ratio_Delta_star_d1=c(1), l=c(0.7962, 2.5204), 
#u=c(2.7625, 2.5204), theta=c(0,0.2,0), pow=0.9, ordering=FALSE, increasing_theta=FALSE, 
#seed=140691, n_trials=10000000, rule=2)
}
