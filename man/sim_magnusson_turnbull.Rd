\name{sim_magnusson_turnbull}
\alias{sim_magnusson_turnbull}
\title{Simulations of trials with GSED}
\description{
\code{sim_magnusson_turnbull} is used to simulate clincal trials with GSED for different type of outcome (survival, binary, continuous).
}
\usage{
sim_magnusson_turnbull(K_stages, N_subsets, f, l, u, ratio_Delta_star_d1, type_outcome, 
param_outcome=NA, n_max=NA, incl_rate=NA, mean_cur_c=NA, HR=NA, nb_required=NA, 
nmax_wait=+Inf, ordering, increasing_theta=FALSE, nsim=1000, seed=42)
}
\arguments{  
  \item{K_stages}{Integer indicating the number of stages in the design.}
  \item{N_subsets}{Integer representing the number of possible subgroups.}
  \item{f}{Vector containing the prevalence rates of each subgroup. Must be of length \code{N_subsets}.}
  \item{l}{Vector containing the lower boundaries for stagewise decisions. Must be of length \code{K_stages}.}
  \item{u}{Vector containing the upper boundaries for stagewise decisions. Must be of length \code{K_stages}.}
  \item{ratio_Delta_star_d1}{Vector containing the ratio between the (observed Fisher) information increments at each stage >1 with the (observed Fisher) information at stage 1. Must be of length \code{K_stages}-1.}
  \item{type_outcome}{A string containing the type of outcome, either "survival", "binary", or "continuous".}
  \item{param_outcome}{Must be supplied only if \code{type_outcome} is equal to "binary" or "continuous". The parameters supplied for the binary outcome must be a list of one element containing a matrix of size 2x\code{N_subsets}. The parameters supplied for the continuous outcome must be a list of two elements containing two matrices of size 2x\code{N_subsets}. The matrices should contain probabilities of response, or the means and variances respectively, for in row control or treatment, and in column the subgroup number.}
  \item{n_max}{Integer representing the maximum number of patients to enroll in a trial. Must be supplied only if \code{type_outcome} is equal to "binary" or "continuous", will be ignored otherwise.}
  \item{incl_rate}{Number representing the inclusion rate. Must be supplied only if \code{type_outcome} is equal to "survival", will be ignored otherwise.}
  \item{mean_cur_c}{Number representing the median survival for the control group. Must be supplied only if \code{type_outcome} is equal to "survival", will be ignored otherwise.}
  \item{HR}{Vector containing the expected hazard ratios for each subgroup. Must be of length \code{N_subsets}. Must be supplied only if \code{type_outcome} is equal to "survival", will be ignored otherwise.}
  \item{nb_required}{Integer indicating the maximum number of events required. Must be supplied only if \code{type_outcome} is equal to "survival", will be ignored otherwise.}
  \item{nmax_wait}{For \code{type_outcome} equal to "survival" only, will be ignored otherwise. If specified, maximum number of patients to include in the trial, the inclusions will be stopped when this
  number is achieved and trial will pursue until the number of events required is achieved. Must be 
  superior to \code{nb_required}. Default value is +Inf.}
  \item{ordering}{Boolean indicating if the subgroups (theta) are ordered.}
  \item{increasing_theta}{Boolean indicating if greater values of theta parameters represent better treatment effects. The default value is set at FALSE.}
  \item{nsim}{Integer indicating the number of trials to simulate. The default value is set at 1000.}
  \item{seed}{Interger representing the seed. The default value is set at 42.}
}
\value{
A list is returned composed of: 
\item{prob_rejec}{Percentage of simulated trials (estimated probability) to reject any subgroup.}
\item{prob_accep}{Percentage of simulated trials (estimated probability) to accept the null hypothesis, that is there is no treatment effect in any subgroup.}
\item{list_keep}{A list of the different subgroups that were selected across all simulated trials.}
\item{pct_keep}{Percentage of selection of each subgroup of \code{list_keep} across all simulated trials.}
\item{rejec_stage}{Vector of percentage of simulated trials (estimated probability) to reject any subgroup at each stage.}
\item{accep_stage}{Vector of percentage of simulated trials (estimated probability) to accept the null hypothesis (that is there is no treatment effect in any subgroup) at each stage.}
\item{mean_pat}{Mean number of patients included across all simulated trials.}
\item{mean_duration}{If \code{type_outcome} is equal to "survival", the trial mean duration across all simulated trials is also returned.}
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
#For testing purpose only, larger number of simulations required (see in comments below)
sim_magnusson_turnbull(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), l=c(0.7962, 2.5204), 
u=c(2.7625, 2.5204), ratio_Delta_star_d1=c(1), type_outcome="binary", param_outcome=
list(matrix(c(0.4,0.4,0.4,0.6,0.6,0.6),nrow=2,ncol=3,byrow=TRUE)), n_max=1496, 
ordering=FALSE, nsim=2, seed=42)

#sim_magnusson_turnbull(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), l=c(0.7962, 2.5204), 
#u=c(2.7625, 2.5204), ratio_Delta_star_d1=c(1), type_outcome="binary", param_outcome=
#list(matrix(c(0.4,0.4,0.4,0.6,0.6,0.6),nrow=2,ncol=3,byrow=TRUE)), n_max=1496, 
#ordering=FALSE, nsim=1000, seed=42)

#sim_magnusson_turnbull(K_stages=2, N_subsets=3, f=c(0.6,0.2,0.2), l=c(0.7962, 2.5204), 
#u=c(2.7625, 2.5204), ratio_Delta_star_d1=c(1), type_outcome="binary", param_outcome=
#list(matrix(c(0.5,0.5,0.5,0.5,0.5,0.5),nrow=2,ncol=3,byrow=TRUE)), n_max=1496, 
#ordering=FALSE, nsim=1000, seed=42)

#sim_magnusson_turnbull(K_stages=2, N_subsets=4, f=c(0.25,0.25,0.25,0.25), l=c(0.98,2.35), 
#u=c(2.59,2.35), ratio_Delta_star_d1=c(1), type_outcome="survival", incl_rate=1/28, 
#mean_cur_c=7/log(2), HR=c(0.8,0.8,0.8,0.8), nb_required=1030, ordering=TRUE, 
#increasing_theta=FALSE, nsim=1000, seed=42)
}
