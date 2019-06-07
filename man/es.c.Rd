\name{es.c}
\alias{es.c}
\title{es.c}
\usage{
es.c(data1,data2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
}
\description{
Calculates the effect size of the difference, its variance and confidence interval from a dataset and a constant.

Default function returns unbiased c (unbiased=TRUE).
When unbiased=FALSE, it returns biased c.
One of the sample size needs to be 1 and it is treated as a constant.
The confidence interval is given by binary search. Therefore, the confidence interval has rounding error, and alpha with many digits will need very long computation time.
In some situations, the confidence interval may be inaccurate, and the warning message appears because of the limitation of the function pt(), which is used for the calculation.
}
\arguments{
\item{data1}{Dataset 1 as a numeric vector.}
\item{data2}{Dataset 2 as a numeric vector.}
\item{alpha}{Confidence level, or type I error rate for the confidence interval. The default value gives 95\% CI.}
\item{unbiased}{When true (default), the bias corrected value is returned. When false, the value based on the whole population sample is returned.}
\item{vector_out}{Whether output is in vector or not. When in vector, it is c(effect size, the variance, lower CI bound, higher CI bound).}
}
\examples{
data1<-c(1,2,3,4)
data2<-c(3)
es.c(data1,data2)
}
\references{
Aoki S., Ito M. & Shimada M. (2019) Effect sizes of the differences between means without assuming the variance equality and between a mean and a constant. arXiv:1901.09581.

Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd edition. New York: Academic Press.
}
