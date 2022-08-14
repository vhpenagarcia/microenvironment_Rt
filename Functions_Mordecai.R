#### Fit existing model developed by Mordecai
# a: Briere
# EFD: Briere
# pEA: Quadratic
# MDR: Briere
# lf: Quadratic
# b: Briere
# c: Briere
# PDR: Briere

##Create a matrix with constant parameters
T0 <- c(13.35,14.58,13.56,11.36,9.16,17.05,12.22,10.68)
Tm <- c(40.08,34.61,38.29,39.17,37.73,35.83,37.46,45.90)
c <- c(2.02e-04,8.56e-03,-5.99e-03,7.86e-05,-1.48e-01,8.49e-04,4.91e-04,6.65e-05)
col_names <- c("T0","Tm","c")
row_names <- c("a","EFD","pEA","MDR","lf","b","c","PDR")
par <- matrix(c(T0,Tm,c),nrow=8,ncol=3,dimname=list(row_names,col_names))

## Briere function
briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b
}

## Quadratic function
quad.2<-function(t, qd, Tm, T0){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<-qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
  }
  b
}

## R0 equation
R0 <- function(a,EFD,pEA,MDR,lf,b,c,PDR){
  mu=1/lf
  ((a^2*b*c*exp(-mu/PDR)*EFD*pEA*MDR)/(mu^3))^0.5
}


## Function to re-scale
rescale <- function(x,min_allowed, max_allowed, min, max){
  if(is.null(min)){min = min(x)} else(min = min)
  if(is.null(max)){max = max(x)} else(max = max)
  ((max_allowed-min_allowed)*(x-min))/((max-min)) + min_allowed
}

## Obtain cases from R0
temp_cases_mor <- function(Temp, sd_T, cases_initial){
  b = c()
  Temp <- append(Temp,25,0)
  sd_T <- append(sd_T,1,0)
  b[1] <- cases_initial
  for(i in 2:length(Temp)){
    temps <- rnorm(1, Temp[i], sd_T[i])
    a_t <- briere(temps,par[1,3],par[1,2],par[1,1])
    b_t <- briere(temps,par[6,3],par[6,2],par[6,1])
    MDR_t <- briere(temps,par[4,3],par[4,2],par[4,1])
    PDR_t <- briere(temps,par[8,3],par[8,2],par[8,1])
    EFD_t <- briere(temps,par[2,3],par[2,2],par[2,1])
    pEA_t <- quad.2(temps,par[3,3],par[3,2],par[3,1])
    lf_t <- quad.2(temps,par[5,3],par[5,2],par[5,1])
    c_t <- briere(temps,par[7,3],par[7,2],par[7,1])
    R0_t <- R0(a=a_t, EFD=EFD_t, pEA=pEA_t, MDR=MDR_t, lf=lf_t, b=b_t, c=c_t, PDR=PDR_t)
    R0_rsc <- rescale(R0_t, 0.01, 1.38, 0, 24.47)
    b[i] <- R0_rsc*b[i-1]
  }
  return(b[-1])
}

