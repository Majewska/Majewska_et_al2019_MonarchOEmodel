## 
## Code to accompany: Multiple transmission routes sustain high prevalence of a virulent parasite in a butterfly host 
## (Majewska, Sims, Schneider, Altizer, Hall 2019, Proceedings B;  DOI: 10.1098/rspb.2019.1630) 
## Code was written by RH and AAM; Latin Hypercube code modified from code provided by DJBecker

####last updated 7.26.19

# Part 2 Latin Hypercube Sampling 

###clear all workspace
rm(list=ls()) 
graphics.off()

## packages
library(diagram)
library(png)
library(deSolve)
library(lhs)
library(sensitivity)
library(betareg)
library(visreg)

### OE-Monarch Model
Monarch_dis = function(t,y,p){
  
  if (t > tau1) Lag1 =lagvalue(t-tau1) else Lag1 = ystart # time lag due to egg stage
  if (t > tau2) Lag2 =lagvalue(t-tau2) else Lag2 = ystart # time lag due to pupa stage
  
  # Variables
  S_L = y[1] # susceptible larvae
  I_V = y[2] # infected larvae via vertical transmission
  I_E = y[3] # infected larvae via environmental transmission
  I_H = y[4] # infected larvae via adult transfer transmission
  S_A = y[5] # susceptible adults
  I_A = y[6] # infected adults
  C_A = y[7] # contaminated adults
  M = y[8]   # milkweed
  W = y[9]   # contaminated milkweed
  
  
  #define ODEs
  with(as.list(params),
       {
         # Total number of larvae:
         Ltot = sum(c(S_L, I_V, I_E, I_H)) 
         
         # Susceptible larvae:
         dS_L = b_s*p_se*Lag1[5] + b_s*p_se*(1-p_h)*Lag1[7] + b_i*p_se*(1-p_v)*Lag1[6] - g*S_L - (mu_d+mu_1*Ltot/M)*S_L - (c*W*S_L/M)
         
         # Infected larvae via vertical transmission
         dI_V = b_i*p_se*p_v*Lag1[6] - g*I_V - (mu_d+mu_1*Ltot/M)*I_V
         
         # Infected larvae via environmental transmission
         dI_E = c*W*S_L/M - (mu_d+mu_1*(Ltot)/M)*I_E -g*I_E
         
         # Infected larvae via adult transfer transmission
         dI_H = b_s*p_se*p_h*Lag1[7] - (mu_d+mu_1*(Ltot)/M)*I_H -g*I_H
         
         # Susceptible adults
         dS_A = g*Lag2[1]*rho-mu_S*S_A- delta*(S_A)*(I_A/(S_A+I_A+C_A)) + mu_C*C_A
         
         # Infected adults
         dI_A = theta*g*(Lag2[2]+Lag2[3]+Lag2[4])*rho-mu_I*I_A
         
         # Contaminated adults due to adult transfer
         dC_A <- delta*(S_A)*(I_A/(S_A+I_A+C_A)) - mu_S*C_A - mu_C*C_A
         
         # Milkweed growth
         dM = r*M*(1-M/K) - c*(Ltot)
         
         # Contaminated milkweed
         dW = lambda*(1-W/M)*(I_A) - mu_W*W - c*(Ltot)*W/M
         
         dy=c(dS_L,dI_V,dI_E,dI_H,dS_A,dI_A,dC_A,dM,dW)   
         list(dy)   
         
       })
}

# Define model parameters

b_s = 9.208333          # Susceptible host fecundity rate (eggs/adult/day)
b_i = 6.5195            # Infected host fecundity rate (eggs/adult/day)
tau1=3                  # Time in egg stage
p_se = 0.5              # probability of egg surviving
g = 1/9                 # Larval development rate 
rho = 0.76              # Probability of pupa surviving to adult
mu_S = 1/24             # Mortality of uninfected adult
mu_I = 1/20             # Mortality rate of infected adult
tau2=7                  # Time in pupa stage
theta = 0.72            # Probability of infected adult eclosing and mating
p_v = 0.9               # Probability of vertical transmission
delta = 1/2             # Daily adult mating probability  
r = 0.032               # Milkweed growth rate
K= 40000                # End of season number of milkweed leaves (assuming 200 leaves per plant)
c = 35/9                # Larval consumption rate of milkweed
mu_W = 0.0125           # Spore decay rate on milkweed 



#######################
##### Starting values
#######################

S_L=0
I_V=0
I_E=0
I_H=0
S_A=18
I_A=2
C_A=0
M=7000
W=0

ystart = c(S_L,I_V,I_E,I_H,S_A,I_A,C_A,M,W) ## the order here matters 

## time
tmax=150
times=seq(0,tmax,by=1)


## ph Probability of contaminated adult infecting their offspring
phmin=0.1
phmax=1

## 1/ mu_C Rate of spore loss for contaminated adults
dmucmin=1
dmucmax=15

## lambda Milkweed leaf visitation rate by infected adults(leaves/day)
lambdamin=1
lambdamax=200

## mu_d Density-independent per capita larval mortality rate 
mu_dmin=0.120
mu_dmax=0.234

## larvdens Larval density
larvdensmin=0.25
larvdensmax=3

## bs Susceptible host fecundity rate (eggs/adult/day)
bsmin=5.7/2
bsmax=15/2

# 1/ mu_W Spore decay rate on milkweed 
dmuWmin=1
dmuWmax=80


## set number of reps
samples=500

## set unknown parameters
upars=7

## lhs sampling
set.seed(5)
lhssample=randomLHS(samples,upars)

## uniform distributions of unknown parameters
bss=(bsmax-bsmin)*lhssample[,1]+bsmin; 
p_hs=(phmax-phmin)*lhssample[,2]+phmin;
dmu_Cs=(dmucmax-dmucmin)*lhssample[,3]+dmucmin; 
lambdas=(lambdamax-lambdamin)*lhssample[,4]+lambdamin; 
mu_ds=(mu_dmax-mu_dmin)*lhssample[,5]+mu_dmin; 
larvdens=(larvdensmax-larvdensmin)*lhssample[,6]+larvdensmin; 
dmu_Ws=(dmuWmax-dmuWmin)*lhssample[,7]+dmuWmin; 


### Calculate mu_1 (density-dependent per capita larval mortality rate)
mu1=function(l=larvden,b_s=b_s,g=g,rho=rho,mu_S=mu_S,mu0=mu_d){
p=(200/l)*((b_s*g*rho/mu_S) - mu0 - g)
return(p)
}

#Calculate number of days a spore-contaminated adult retains spores 
muC=function(d=dmu_C){
  p=1/d
  return(p)
}

#Calculate number of days a spore-contaminated milkweed leaf remains infectious
muW=function(w=dmu_W){
  p = 1/w
  return(p)
}

## set empty vectors
ps=rep(NA,samples)
mu_1s=rep(NA,samples)
mu_Cs=rep(NA,samples)
mu_Ws=rep(NA,samples)

## loop
for(nsample in 1:samples) ## initiate loop
{
  
  ## start loop
  print(sprintf('Starting Simulation %d of %d',nsample,samples));
  
  ## values for lhs parameters
  b_s = bss[nsample]; print(sprintf('b_s = %f',b_s));
  p_h = p_hs[nsample]; print(sprintf('p_h = %f',p_h));
  dmu_C = dmu_Cs[nsample]; print(sprintf('dmu_C = %f',dmu_C));
  dmu_W = dmu_Ws[nsample]; print(sprintf('dmu_W = %f',dmu_W));
  lambda = lambdas[nsample]; print(sprintf('lambda = %f',lambda));
  mu_d=mu_ds[nsample]; print(sprintf('mu_d = %f',mu_d));
  larvden=larvdens[nsample]; print(sprintf('larvden = %f',larvden));
  
  mu_1s[nsample]=mu1(l=larvden,b_s=b_s,g=g,rho=rho,mu_S=mu_S,mu0=mu_d)
  mu_1=mu_1s[nsample]; print(sprintf('mu_1 = %f',mu_1))
  
  mu_Cs[nsample]=muC(d=dmu_C)
  mu_C=mu_Cs[nsample]; print(sprintf('mu_C = %f',mu_C))
  
  mu_Ws[nsample]=muW(w=dmu_W)
  mu_W=mu_Ws[nsample]; print(sprintf('mu_W = %f',mu_W))
  
 
  ## set all parameters
  params = c(r=r, K=K, b_s=b_s,b_i=b_i,mu_d=mu_d,mu_1=mu_1,g=g,c=c,theta=theta,lambda=lambda,p_v=p_v,mu_I=mu_I,mu_S=mu_S,mu_W=mu_W,mu_C=mu_C,delta=delta,p_h=p_h,tau1=tau1,tau2=tau2)
  
  ## run ode
  out = as.data.frame(dede(func=Monarch_dis,y=ystart,t=times, p=params))
  names(out)=c("times","S_L","I_V","I_E","I_H","S_A","I_A","C_A","M","W") ## make names easy to call
  
  out$noa=with(out,S_A+I_A+C_A) ## calculate number of adults
  out$prev=with(out,(I_A)/noa) ## calculate total proportion of infected adults

  ## clean values
  out$noa=ifelse(out$noa<0,0,out$noa) ## make zero if less than zero (weirdness)
  out$prev=ifelse(out$prev<0,0,out$prev) ## make zero if less than zero (weirdness)
  
  ps[nsample]=tail(out$prev,1)

}

## make parameter data frame
lhsdata=data.frame(bss,mu_ds,larvdens,p_hs,dmu_Cs,lambdas,dmu_Ws,ps)

names(lhsdata)=c("bs","mu0","larvdens","ph","dmuC","lambd","dmuW","prevstar")

inputs=lhsdata

inputs$prevstar=NULL


## Partial Rank Correlation Coefficient (PRCC) of prevalance (here 'ps')
set.seed(1)
cordat=pcc(X=inputs,y=ps,rank=T,conf=0.95,nboot=100)
cordat2=round(cordat$PRCC,2)
cordat2$par=rownames(cordat2)
pstarcor=cordat2

#Labels
labs=c(expression(paste(b[s])),
       expression(paste(mu[0])),
       expression(paste(L)), 
       expression(paste(p[h])),
       expression(paste(T[C])),
       expression(paste(lambda)),
      expression(paste(T[W])))

## Plot space

lhsdata$pstar2=((lhsdata$prevstar*(nrow(lhsdata)-1))+0.5)/nrow(lhsdata)

#Figure S3a
#png("FigureS3a.png",width=8,height=3,units="in",res=150)
par(mfrow=c(1,ncol(inputs)),mar=c(1.5,1.5,1,1),oma=c(1,4,1,1))
pmod=betareg(pstar2~bs+ph+dmuC+larvdens+lambd+mu0+dmuW,data=lhsdata)

## first plot for pstar
plot(lhsdata$bs,lhsdata$pstar2,las=1,xlab="",ylab="",pch=21,bg="gray",col="gray")
mtext(labs[1],side=3,line=0.3)

## add visreg fit
fit=visreg(pmod,"bs",plot=F,scale="response")
lines(visregFit~bs,data=fit$fit,lwd=2)

## add mtext for P*
text(-1,.5, labels = 'Late season infection prevalence', xpd = NA, srt = 90, cex = 1.25)

## loop through
base=inputs[-1]
labs2=labs[-1]
scor=pstarcor[-1,]
for(i in 1:ncol(base)){
  
  ## sub data
  set=data.frame(base[i],lhsdata["pstar2"])
  
  ## plot and label
  plot(set,las=1,xlab="",ylab="",yaxt="n",pch=21,bg="gray",col="gray")
  mtext(labs2[i],side=3,line=0.3)
  
  ## visreg fit
  fit=visreg(pmod,names(set)[1],plot=F,scale="response")
  test=fit$fit
  test=test[c(names(set)[1],"visregFit")]
  lines(test,lwd=2)
  
}
#dev.off()

graphics.off()

#Figure S3b
## Correlations only
#png("FigureS3b.png",width=8,height=5,units="in",res=300)
plot(0,0,type="n",las=1,ylim=c(-1,1),xlim=c(0.5,nrow(pstarcor)+0.5),
     xlab="",xaxt="n",cex.lab=1.25,
     ylab=expression(paste("PRCC with late season infection prevalence")))
abline(v=1:nrow(pstarcor),lty=3)
abline(h=0,lty=2,lwd=2)
segments(x0=1:nrow(pstarcor),x1=1:nrow(pstarcor),
         y0=pstarcor$`min. c.i.`,y1=pstarcor$`max. c.i.`,lwd=4)
points(1:nrow(pstarcor),pstarcor$original,pch=21,bg="grey",cex=2,lwd=2)

## axis
axis(1,at=1:nrow(pstarcor),labels=labs,cex.axis=1.55)

#dev.off()

