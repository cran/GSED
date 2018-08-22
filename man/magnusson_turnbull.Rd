\name{magnusson_turnbull}
\alias{magnusson_turnbull}
\title{Application of GSED on data}
\description{
\code{magnusson_turnbull} is used apply GSED design, selection or evaluation at each stage, on data.
}
\usage{
magnusson_turnbull(stage_cur, keep=NA, N_subsets, Y, I, l, u, ordering, 
increasing_theta=FALSE)
}
\arguments{                
  \item{stage_cur}{Integer representing the current stage. 0 represents selection at stage 1, 1 represents evaluation at stage 1, while k (>1) represents evaluation at stage k.}
  \item{keep}{Vector of indices of selected subgroups if selection at stage 1 is already performed. Values must be between 1 and \code{N_subsets}. By default filled with NA if the function is run for selection step.}
  \item{N_subsets}{Integer representing the number of possible subgroups.}
  \item{Y}{Efficient score test statistics. For \code{stage_cur}>0 (evaluation at stage 1 or k (k>1),), value representing the efficient score test statistic for all (pooled) selected subgroup. For \code{stage_cur}=0 (selection at stage 1), vector representing the efficient score test statistic for each subgroup.}
  \item{I}{Observed Fisher information. For \code{stage_cur}>0 (evaluation at stage 1 or k (k>1),), value representing the observed Fisher information for all (pooled) selected subgroup. For \code{stage_cur}=0 (selection at stage 1), vector representing the observed Fisher information for each subgroup.}
  \item{l}{Vector containing the lower boundaries for stagewise decisions.}
  \item{u}{Vector containing the upper boundaries for stagewise decisions.}
  \item{ordering}{Boolean indicating if the subgroups (theta) are ordered.}
  \item{increasing_theta}{Boolean indicating if greater values of theta parameters represent better treatment effects. The default value is set at FALSE.}
}
\value{
An list is returned, consisting of: 
\item{Rejection}{Interger with value 1 if the decision is to reject the null hypothesis, 0 otherwise.}
\item{Acceptation}{Interger with value 1 if the decision is to accept the null hypothesis, 0 otherwise.}
\item{Keep}{Vector of indices of selected subgroups (between 1 and \code{N_subsets}).}
}
\references{
Baldur P. Magnusson and Bruce W. Turnbull. Group sequential enrichment design incorporating subgroup selection. Statistics in Medicine, 2013. <doi:10.1002/sim.5738>
}
\author{Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\examples{
magnusson_turnbull(stage_cur=0, keep=NA, N_subsets=3, Y=c(-10.71,12.84,19.06), 
I=c(480,144,176), l=c(0.7962,2.5204), u=c(2.7625,2.5204), ordering=FALSE)

magnusson_turnbull(stage_cur=2, keep=c(2,3), N_subsets=3, Y=135.57, 
I=1120, l=c(0.7962,2.5204), u=c(2.7625,2.5204), ordering=FALSE)
}