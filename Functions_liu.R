#### Equations for Liu-Helmersson 2014 model

### Just temperatures
## Biting rate a
a <- function(T){
  0.0043*T + 0.0943
}


##Probability of infection from humans to vectors per bite: bm
bm <- function(T){
  b=c()
  for(i in 1:length(T)){
    if(T[i]<12.4) {b[i]<-0}
    else if(T[i]>12.4 && T[i]<26.1) {b[i]<-0.0729*T[i]-0.9037}
    else if(T[i]>26.1) {b[i]<-1}
  }
  b
}

## The probability of transmission from vector to human per bite: bh
bh <- function(T){
  b=c()
  for(i in 1:length(T)){
    if(T[i]<12.286) {b[i] <- 0}
    else if(T[i]>12.286 && T[i]<32.51) # Modified upper limit from 32.461 to 32.51
      {b[i] <- 0.001044*T[i]*(T[i]-12.286)*sqrt(32.51-T[i])}
    else if(T[i]>32.51) {b[i] <- 0}
  }
  b
}


## Extrinsic incubation Period: n
n <- function(T){
  4 + exp(5.15-0.123*T)
}


## Mortality rate: um
um <- function(T){
  mort <- 0.8692-0.1590*T+0.01116*(T^2)-3.408e-4*(T^3)+3.809e-6*(T^4)
  um_cond <- ifelse(mort >= 1, 1, mort)
  return(um_cond)
}


## Vectorial capacity
rVC <- function(a,bh,bm,um,n){
  rvc <- (a^2 * bh * bm * exp(-um*n))/um
  rvc2 <- ifelse(rvc == 0, 0.01, rvc)
  return(rvc2)
}



## Function to return cases directly from temp
temp_cases_liu <- function(Temp, sd_T, cases_initial){
  b = c()
  Temp <- append(Temp,25,0)
  sd_T <- append(sd_T,1,0)
  b[1] <- cases_initial
  for(i in 2:length(Temp)){
    temps <- rnorm(1, Temp[i], sd_T[i])
    a_t <- a(temps)
    bh_t <- bh(temps)
    bm_t <- bm(temps)
    um_t <- um(temps)
    n_t <- n(temps)
    rvc_t <- rVC(a=a_t,bh=bh_t,bm=bm_t,um=um_t,n=n_t)
    b[i] <- rvc_t*b[i-1]
  }
  return(b[-1])
}
