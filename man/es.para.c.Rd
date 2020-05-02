\name{es.para.c}
\alias{es.para.c}
\title{es.para.c}
\usage{
es.para.c(mean1,mean2,var1,var2,n1,n2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
}
\description{
Calculates the effect size of the difference, its variance and confidence interval from statistics of one dataset and a constant.
Default function returns unbiased c (unbiased=TRUE).
When unbiased=FALSE, it returns biased c.
One of the sample size needs to be 1 and it is treated as a constant. The confidence interval is given by binary search. Therefore, the CI has rounding error, and alpha with many digits will need very long computation time.
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
es.para.c(2.5,3,5/3,0,4,1)
}
\references{
Aoki S. (2020) Effect sizes of the differences between means without assuming the variance equality and between a mean and a constant. Heliyon 6(1): e03306.

Aoki S., Ito M. & Shimada M. (2019) Effect sizes of the differences between means without assuming the variance equality and between a mean and a constant. arXiv:1901.09581.

Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd edition. New York: Academic Press.
}
