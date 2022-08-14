#### Function by Caminade 2017

## Extrinsic incubation period
a1 <- function(T){
  0.0043*T+0.0943
}

a2 <- function(T){
  a1 <- 0.0043*T + 0.0943
  a2 <- 0.5*a1
  return(a2)
}

##Vector preferences
phi1 <- 1
phi2 <- 0.5

## Transmission probability -vector to host
b1 <- 0.5
b2 <- 0.5

## Transmission probability -host to vector
beta1 <- 0.1
beta2 <- 0.033

##Mortality rates
mu1 <- function(T){
  b=c()
  for(i in 1:length(T)){
    if(T[i]<22) {b[i] <- 1/(1.22 + exp(-3.05 + 0.72*T[i]))+0.196}
    else if(T[i]>=22) {b[i] <- 1/(1.14 + exp(51.4 - 1.3*T[i]))+0.192}
  }
  b
}

mu2 <- function(T){
  b=c()
  for(i in 1:length(T)){
    if(T[i]<15) {b[i] <- 1/(1.1 + exp(-4.04 + 0.576*T[i]))+0.12}
    else if(T[i]>=15 && T[i]<26.3) {b[i] <- 0.000339*T[i]^2 - 0.0189*T[i] + 0.336}
    else if(T[i]>=26.3) {b[i] <- 1/(1.065 + exp(32.2 - 0.92*T[i]))+0.0747}
  }
  b
}

## EIP (days)= 1/v, v=1/EIP
v1 <- function(T){
  1/(4 + exp(5.15-0.123*T))
}

v2 <- function(T){
  1/(1.03*(4 + exp(5.15-0.123*T)))
}


##Recovery rate (per day)
r <- 1/7

### R for each Aedes species

R <- function(b,beta,a,mu,v,phi,m,r){
  R_1 <- (b*beta*a^2)/mu
  R_2 <- v/(v+mu)
  R_3 <- (phi^2*m)/r
  return(R_1*R_2*R_3)
}

###R0
R0 <- function(R11,R22){
  sqrt(R11+R22)
}

### Rescale R0
rescale <- function(x,min_allowed, max_allowed, min, max){
  if(is.null(min)){min = min(x)} else(min = min)
  if(is.null(min)){max = max(x)} else(max = max)
  ((max_allowed-min_allowed)*(x-min))/((max-min)) + min_allowed
}

### Obtain cases from temperature
temp_cases_cam <- function(Temp, sd_T, cases_initial, m1, m2){
  b = c()
  Temp <- append(Temp,25,0)
  sd_T <- append(sd_T,1,0)
  b[1] <- cases_initial
  for(i in 2:length(Temp)){
    temps <- rnorm(1, Temp[i], sd_T[i])
    a1_t <- a1(temps)
    a2_t <- a2(temps)
    mu1_t <- mu1(temps)
    mu2_t <- mu2(temps)
    v1_t <- v1(temps)
    v2_t <- v2(temps)
    R01_t <- R(b=b1,beta=beta1,a=a1_t,mu=mu1_t,v=v1_t, phi=phi1,m=m1,r=r)
    R02_t <- R(b=b2,beta=beta2,a=a2_t,mu=mu2_t,v=v2_t, phi=phi2,m=m2,r=r)
    R0_t <- R0(R11=R01_t, R22=R02_t)
    R0_rsc <- rescale(R0_t, 0.01, 1.38, 1, 6.4)
    b[i] <- R0_rsc*b[i-1]
  }
  return(b[-1])
}

