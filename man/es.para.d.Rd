\name{es.para.d}
\alias{es.para.d}
\title{es.para.d}
\usage{
es.para.d(mean1,mean2,var1,var2,n1,n2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
}
\description{
Calculates the effect size of the difference, its variance and confidence interval from statistics of two datasets under the variance equality.
Default function returns Hedges' d (unbiased=TRUE).
When unbiased=FALSE, it returns Cohen's d. The confidence interval is given by binary search. Therefore, the confidence interval has rounding error, and alpha with many digits will need very long computation time.
In some situations, the confidence interval may be inaccurate, and the warning message appears because of the limitation of the function pt(), which is used for the calculation.
}
\arguments{
\item{mean1}{Mean of dataset 1.}
\item{mean2}{Mean of dataset 2.}
\item{var1}{Unbiased (divided by n1-1) variance of dataset 1.}
\item{var2}{Unbiased (divided by n2-1) variance of dataset 2.}
\item{n1}{Sample size of dataset 1.}
\item{n2}{Sample size of dataset 2.}
\item{alpha}{Confidence level, or type I error rate for the confidence interval. The default value gives 95\% CI.}
\item{unbiased}{When true (default), the bias corrected value is returned. When false, the value based on the whole population sample is returned.}
\item{vector_out}{Whether output is in vector or not. When in vector, it is c(effect size, the variance, lower CI bound, higher CI bound).}
}
\examples{
es.para.d(2,3,1,2,4,3)
}
\references{
Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd edition. New York: Academic Press.

Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators.
Journal of Educational Statistics 6, 107--128.

Hedges, L. V. and Olkin I. (1985). Statistical Methods for Meta-analysis. Orlando: Academic Press.
}
