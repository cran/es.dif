es.d<-function(data1,data2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-inputCheck(data1,data2,alpha,unbiased,vector_out)
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  data1<-data1[!is.na(data1)]
  data2<-data2[!is.na(data2)]
  n1<-length(data1)
  n2<-length(data2)
  mean1<-mean(data1)
  mean2<-mean(data2)
# compute the pooled standard deviation
  if(n1!=1)
  {
    var1<-stats::var(data1)
  }
  else
  {
    var1<-0
  }
  if(n2!=1)
  {
    var2<-stats::var(data2)
  }
  else
  {
    var2<-0
  }
  return(es.para.d(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out))
}

es.para.d<-function(mean1,mean2,var1,var2,n1,n2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-input.para.Check(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out)
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  if(n1+n2<=3)
  {
    stop("Error:The sum of the sample size must be 4 or larger.")
  }
  # compute the pooled standard deviation
  ss1<-var1*(n1-1)
  ss2<-var2*(n2-1)
  df<-n1+n2-2
  div<-sqrt((ss1+ss2)/df)
  # compute the correction factor J
  J<-calcJ(df)
  # compute d
  if(unbiased)
  {
    d<-(mean1-mean2)/div*J
  }
  else
  {
    d<-(mean1-mean2)/div
  }
  # compute the variance
  enye<-n1*n2/(n1+n2)
  if(unbiased)
  {
    para.d<-d
  }
  else
  {
    para.d<-d*J
  }
  var.d<-df/(df-2)*(1/enye+para.d^2)-para.d^2/J^2
  if(unbiased)
  {
    var.d<-var.d*J^2
  }
  #compute CI
  t<-(mean1-mean2)/div*sqrt(enye)
  CI<-findCI(alpha,t,df)
  if(unbiased)
  {
    CI<-CI/sqrt(enye)*J
  }
  else
  {
    CI<-CI/sqrt(enye)
  }
  #output
  if(vector_out)
  {
    return(c(d,var.d,CI))
  }
  else
  {
    result <- matrix(1, nrow=3,ncol=2)
    if(unbiased)
    {
      result[1,1]<-"Hedges' d:"
    }
    else
    {
      result[1,1]<-"Cohen's d:"
    }
    result[1,2]<-d
    result[2,1]<-"variance:"
    result[2,2]<-var.d
    result[3,1]<-"CI:"
    result[3,2]<-paste("[",as.character(CI[1]),",",as.character(CI[2]),"]")
    return(result)
  }
}

es.e<-function(data1,data2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-inputCheck(data1,data2,alpha,unbiased,vector_out)
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  data1<-data1[!is.na(data1)]
  data2<-data2[!is.na(data2)]
  n1<-length(data1)
  n2<-length(data2)
  mean1<-mean(data1)
  mean2<-mean(data2)
  # compute the pooled standard deviation
  if(n1!=1)
  {
    var1<-stats::var(data1)
  }
  else
  {
    var1<-0
  }
  if(n2!=1)
  {
    var2<-stats::var(data2)
  }
  else
  {
    var2<-0
  }
  return(es.para.e(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out))
}
es.para.e<-function(mean1,mean2,var1,var2,n1,n2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-input.para.Check(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out)
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  if((n1==1)||(n2==1))
  {
    stop("Error: Datum with the sample size 1 cannot be used.")
  }
  df<-(var1/n1+var2/n2)^2/((var1/n1)^2/(n1-1)+(var2/n2)^2/(n2-1))
  div<-sqrt(var1/n1+var2/n2)
  enye<-n1*n2/(n1+n2)
  # compute the correction factor J
  J<-calcJ(df)
  # compute e
  if(unbiased)
  {
    e<-(mean1-mean2)/div*J/sqrt(enye)
  }
  else
  {
    e<-(mean1-mean2)/div/sqrt(enye)
  }
  # compute the variance
  if(unbiased)
  {
    para.e<-e
  }
  else
  {
    para.e<-e*J
  }
  var.e<-df/(df-2)*(1/enye+para.e^2)-para.e^2/J^2
  if(unbiased)
  {
    var.e<-var.e*J^2
  }
  #compute CI
  t<-(mean1-mean2)/div
  CI<-findCI(alpha,t,df)
  if(unbiased)
  {
    CI<-CI/sqrt(enye)*J
  }
  else
  {
    CI<-CI/sqrt(enye)
  }
  #output
  if(vector_out)
  {
    return(c(e,var.e,CI))
  }
  else
  {
    result <- matrix(1, nrow=3,ncol=2)
    if(unbiased)
    {
      result[1,1]<-"Unbiased e:"
    }
    else
    {
      result[1,1]<-"Biased e:"
    }
    result[1,2]<-e
    result[2,1]<-"variance:"
    result[2,2]<-var.e
    result[3,1]<-"CI:"
    result[3,2]<-paste("[",as.character(CI[1]),",",as.character(CI[2]),"]")
    return(result)
  }
}

es.c<-function(data1,data2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-inputCheck(data1,data2,alpha,unbiased,vector_out)
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  data1<-data1[!is.na(data1)]
  data2<-data2[!is.na(data2)]
  n1<-length(data1)
  n2<-length(data2)
  mean1<-mean(data1)
  mean2<-mean(data2)
  # compute the standard deviation
  if(n1!=1)
  {
    var1<-stats::var(data1)
  }
  else
  {
    var1<-0
  }
  if(n2!=1)
  {
    var2<-stats::var(data2)
  }
  else
  {
    var2<-0
  }
  return(es.para.c(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out))
}

es.para.c<-function(mean1,mean2,var1,var2,n1,n2,alpha=0.05,unbiased=TRUE,vector_out=FALSE)
{
  errorCheck<-input.para.Check(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out)
  if(n1+n2<=3)
  {
    stop("Error: the sample size of non-constant sample must be 3 or larger.")
  }
  if(errorCheck!=0)
  {
    return(errorCheck)
  }
  # compute the standard deviation
  if(n1!=1 && n2==1)
  {
    ss<-var1*(n1-1)
    n<-n1
  }
  else if(n1==1 && n2!=1)
  {
    ss<-var2*(n2-1)
    n<-n2
  }
  else if(n1==1 && n2==1)
  {
    return("Error: both of the data size are 1.")
  }
  else if(n1!=1 && n2!=1)
  {
    return("Error: both of the data size are non-1.")
  }
  df<-n-1
  div<-sqrt(ss/df)
  # compute the correction factor J
  J<-calcJ(df)
  # compute c
  if(unbiased)
  {
    c<-(mean1-mean2)/div*J
  }
  else
  {
    c<-(mean1-mean2)/div
  }
  # compute the variance
  if(unbiased)
  {
    para.c<-c
  }
  else
  {
    para.c<-c*J
  }
  var.c<-df/(df-2)*(1/df+para.c^2)-para.c^2/J^2
  if(unbiased)
  {
    var.c<-var.c*J^2
  }
  #compute CI
  t<-(mean1-mean2)/div*sqrt(n-1)
  CI<-findCI(alpha,t,df)
  if(unbiased)
  {
    CI<-CI/sqrt(n-1)*J
  }
  else
  {
    CI<-CI/sqrt(n-1)
  }
  #output
  if(vector_out)
  {
    return(c(c,var.c,CI))
  }
  else
  {
    result <- matrix(1, nrow=3,ncol=2)
    if(unbiased)
    {
      result[1,1]<-"Unbiased c:"
    }
    else
    {
      result[1,1]<-"Biased c:"
    }
    result[1,2]<-c
    result[2,1]<-"variance:"
    result[2,2]<-var.c
    result[3,1]<-"CI:"
    result[3,2]<-paste("[",as.character(CI[1]),",",as.character(CI[2]),"]")
    return(result)
  }
}

##function to calculate the correction coefficient J
#this function is not for users' direct use
calcJ<-function(m)
{
  J<-1
  if((m>1)&&(m<=171))
  {
    J<-gamma(m/2)/sqrt(m/2)/gamma((m-1)/2)
  }
  else if(m>171)
  {
    J<-(4*m-4)/(4*m-1)
  }
  else if(m==1)
  {
    J<-0
  }
  else
  {
    J<-NaN
    stop("Unexpected Error: negative degree of freedom.")
  }
  return(J)
}

##function to search CIs
#this function is not for users' direct use
findCI<-function(alpha,t,df)
{
  #check for infinit t
  if(is.infinite(t))
  {
    return(c(NaN,NaN))
  }
  #find number of digit of alpha/2
  x<-alpha/2
  digit<-0
  while(x %% 1 != 0)
  {
    x<-x*10
    digit<-digit+1
  }
  #extra digit for more precise computation
  digit<-digit+1
  #search limits
  maxAim<-1-alpha/2
  minAim<-alpha/2
  #search max
  testNcp<--5
  nowP<-stats::pt(t,df,testNcp)

  while(maxAim>=nowP)
  {
    testNcp<-testNcp*2
    nowP<-stats::pt(t,df,testNcp)
  }
  minL<-testNcp

  #search min
  testNcp<-5
  while(minAim<=nowP)
  {
    testNcp<-testNcp*2
    nowP<-stats::pt(t,df,testNcp)
  }
  maxL<-testNcp

  #serach CIs

  result<-numeric(2)
  for(i in 1:2)
  {
    nowMaxL<-maxL
    nowMinL<-minL
    nowNcp<-minL+(maxL-minL)/2
    nowP<-stats::pt(t,df,nowNcp)
    nowP<-signif(nowP,digit)
    if(i==1)
    {
      aimP<-1-alpha/2
    }
    else
    {
      aimP<-alpha/2
    }
    while(nowP!=aimP)
    {
      if(nowP<aimP)
      {
        nowMaxL<-nowNcp
      }
      else
      {
        nowMinL<-nowNcp
      }
      nowNcp<-nowMinL+(nowMaxL-nowMinL)/2
      nowP<-stats::pt(t,df,nowNcp)
      nowP<-signif(nowP,digit)
    }
    result[i]<-nowNcp
  }
  return(result)
}

#function to check input
#this function is not for users' direct use
inputCheck<-function(data1,data2,alpha,unbiased,vector_out)
{
  if(is.character(data1))
  {
    return("Error: data1 contains character.")
  }
  if(is.character(data2))
  {
    return("Error: data2 contains character.")
  }
  if(length(data1)==0)
  {
    return("Error: length of data1 is 0.")
  }
  if(length(data2)==0)
  {
    return("Error: length of data2 is 0.")
  }
  if(is.character(alpha))
  {
    return("Error: alpha must be a number.")
  }
  if(alpha<=0||alpha>1)
  {
    return("Error: alpha must be within (0,1].")
  }
  if(!is.logical(unbiased))
  {
    return("Error: unbiased must be TRUE or FALSE.")
  }
  if(!is.logical(vector_out))
  {
    return("Error: vector_out must be TRUE or FALSE.")
  }
  return(0)
}
#function to check input
#this function is not for users' direct use
input.para.Check<-function(mean1,mean2,var1,var2,n1,n2,alpha,unbiased,vector_out)
{
  if(is.character(mean1))
  {
    return("Error: mean1 contains character.")
  }
  if(is.character(mean2))
  {
    return("Error: mean2 contains character.")
  }
  if(is.character(var1))
  {
    return("Error: var1 contains character.")
  }
  if(is.character(var2))
  {
    return("Error: var2 contains character.")
  }
  if(var1<0)
  {
    return("Error: var1 must be 0 or larger.")
  }
  if(var2<0)
  {
    return("Error: var2 must be 0 or larger.")
  }
  if(is.character(n1))
  {
    return("Error: n1 must be a number.")
  }
  if(is.character(n2))
  {
    return("Error: n2 must be a number.")
  }
  if(n1 %% 1!=0||n1<=0)
  {
    return("Error: n1 must be a natural number.")
  }
  if(n2 %% 1!=0||n2<=0)
  {
    return("Error: n2 must be a natural number.")
  }
  if(is.character(alpha))
  {
    return("Error: alpha must be a number.")
  }
  if(alpha<=0||alpha>1)
  {
    return("Error: alpha must be within (0,1].")
  }
  if(!is.logical(unbiased))
  {
    return("Error: unbiased must be TRUE or FALSE.")
  }
  if(!is.logical(vector_out))
  {
    return("Error: vector_out must be TRUE or FALSE.")
  }
  return(0)
}
