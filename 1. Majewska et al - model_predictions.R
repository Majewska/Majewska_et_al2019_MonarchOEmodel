## 
## Code to accompany: Multiple transmission routes sustain high prevalence of a virulent parasite in a butterfly host 
## (Majewska, Sims, Schneider, Altizer, Hall 2019, Proceedings B;  DOI: 10.1098/rspb.2019.1630) 
## Code was written by RH and AAM and contains the model and main text figures
##

####last updated 7.26.19

# Part 1 Model and model predictions 

###clear all workspace
rm(list=ls()) 
graphics.off()

###load necessary packages
require(deSolve)
require(ggplot2)
require(doBy)
require(gridExtra)
require(reshape2)


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
  M = y[8] # milkweed
  W = y[9] # contaminated milkweed
  
  
  
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
larvdens = 1.3          # Larval density
rho = 0.76              # Probability of pupa surviving to adult
mu_S = 1/24             # Mortality of uninfected adult
mu_I = 1/20             # Mortality rate of infected adult
mu_d = -log(0.06)/9     # Density-independent per capita larval mortality rate 
mu_1 = (200/larvdens)*((b_s*p_se*g*rho/mu_S) - mu_d - g)   
                        # Density-dependent per capita larval mortality rate 
tau2=7                  # Time in pupa stage
theta = 0.72            # Probability of infected adult eclosing and mating
p_v = 0.9               # Probability of vertical transmission
p_h =  0.614            # Probability of contaminated adult infecting their offspring
delta = 1/2             # Daily adult mating probability  
mu_C= 1/14              # Rate of spore loss for contaminated adults
r = 0.032               # Milkweed growth rate
K= 40000                # End of season number of milkweed leaves (assuming 200 leaves per plant)
c = 35/9                # Larval consumption rate of milkweed
lambda = 50             # Milkweed leaf visitation rate by infected adults,  (leaves/day)
mu_W = 0.0125           # Spore decay rate on milkweed 


# Define the range of time over which we want to solve the equations
times = seq(0,150,by=1) 

# Initial conditions
ystart = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=18,I_A=2,C_A=0,M=7000,W=0)

# Combine paramameters into single vector
params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)



################################################################
########## Model predictions ###################################
################################################################

out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)

# Time vector
time<-out[,1]

S_L_out = out[,2] #Number of susceptible larvae over time
I_V_out = out[,3] #Number of infected larvae via vertical transmission over time
I_E_out = out[,4] #Number of infected larvae via environmental transmission over time
I_H_out = out[,5] #Number of infected larvae via adult transfer transmission over time
ITot_out = I_V_out+I_H_out+I_E_out #Number of all infected larvae over time
Ltot_out = S_L_out+I_V_out+I_H_out+I_E_out #Number of all larvae over time

SA_out = out[,6] #Number of susceptible adults over time
IA_out = out[,7] #Number of infected adults over time
CA_out = out[,8] #Number of contaminated adults over time
N_adults= SA_out+CA_out+IA_out #Number of all adults over time

M_out<-out[,9] #Number of milkweed leaves over time
W_out<-out[,10] #Number of contaminated milkweed leaves over time

# Model predictions of:
SA_final = tail(out,n=1)[6] #end of season number of susceptible adults
IA_final = tail(out,n=1)[7] #end of season number of infected adults
CA_final = tail(out,n=1)[8] #end of season number of contaminated adults
NA_final = SA_final+CA_final+IA_final #end of season adult popupaltion size 
prev_final=IA_final/NA_final  # end of season prevalence


###############################
# OE Model Prediction Figures #
###############################

## Model predictions of proportion infected adults over time ##

infA<-IA_out/N_adults
infL<-ITot_out/Ltot_out

IAL<-data.frame(cbind(time,infA))
IAL$stage<-"Adult"

IAL2<-data.frame(cbind(time,infL))
IAL2$stage<-"Larva"

IAL2<-renameCol(IAL2,c("infL"), c("Prop"))
IAL<-renameCol(IAL,c("infA"), c("Prop"))
df<-rbind(IAL2,IAL)

#Figure 2c
ggplot() +
geom_line(data=df,aes(x=time, y=Prop,linetype=stage, color = stage), size = 1.5,na.rm=TRUE, stat="identity")+
theme_classic()+
  labs(x="Time (days)", y = "Proportion infected")+
  theme(axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.key.width = unit(1.5,"cm"))+
  scale_y_continuous(expand = c(0,0),limits = c(0, 1.024))+
  scale_linetype_manual(values = c("solid", "solid"))+
  scale_color_manual(values=c('black','gray65'))+  
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

## Model predictions of contaminated adults and milkweed over time ##

ContA<-CA_out/(CA_out + SA_out)
ContM<-W_out/M_out
C1<-data.frame(cbind(time,ContA))
C1$stage<-"Adults"
C2<-data.frame(cbind(time,ContM))
C2$stage<-"Leaves"
C1<-renameCol(C1,c("ContA"), c("Prop"))
C2<-renameCol(C2,c("ContM"), c("Prop"))
dfM<-rbind(C1,C2)

#Figure 2d
ggplot(dfM, aes(x=Month, y=prop_m))+
  geom_line(data=dfM,aes(x=time, y=Prop,linetype=stage, color = stage), size = 1.5,na.rm=TRUE, stat="identity")+
   theme_classic()+
  labs(fill="", x="Time (days)", y = "Proportion contaminated")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 1.024))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_color_manual(values=c('black','gray'))+
  theme(axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.key.width = unit(1.5,"cm"))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

##############################################################
##### DISEASE - FREE Population Size ######
##############################################################

# Initial conditions
ystart2 = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=20,I_A=0,C_A=0,M=7000,W=0)

# Transmission routes turned off 
p_v = 0               
p_h = 0

params2 = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_free = dede(func=Monarch_dis,y=ystart2,times=times, parms=params2)

# Model predictions of:
SA_final_f = tail(out_free,n=1)[6] #end of season number of susceptible adults
IA_final_f = tail(out_free,n=1)[7] #end of season number of infected adults
CA_final_f = tail(out_free,n=1)[8] #end of season number of contaminated adults
NA_final_f = SA_final_f+CA_final_f+IA_final_f #end of season adult popupaltion size 
prev_final_f=IA_final_f/NA_final_f  # end of season prevalence

NA_final_f

##################################################################################
  ####### OE prediction of importance of different transmission routes ##########
##################################################################################

Vert <- I_V_out/ITot_out
Env <- I_E_out/ITot_out
Hort <- I_H_out/ITot_out

R1 <- data.frame(cbind(time,Vert))
R1$Route <- "V  (vertical transmission)"

R2 <- data.frame(cbind(time,Env))
R2$Route <- "E  (environmental transmission)"

R3 <- data.frame(cbind(time,Hort))
R3$Route <- "A  (adult transfer)"

R1<-renameCol(R1,c("Vert"), c("RelImp"))
R2<-renameCol(R2,c("Env"), c("RelImp"))
R3<-renameCol(R3,c("Hort"), c("RelImp"))
dfRI<-rbind(R1,R2)
dfRI<-rbind(dfRI,R3)

dfRI$Route<-as.factor(dfRI$Route)

dfRI$Route<- factor(dfRI$Route, levels = rev(levels(dfRI$Route)))
  
#Figure 3a
ggplot()+
  geom_line(data=dfRI,aes(x=time, y=RelImp,linetype=Route), size = 1.5,na.rm=TRUE, stat="identity")+
  theme_classic()+
  labs(fill="", x="Time (days)", y = "Relative importance")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))+
  scale_color_manual(values=c('black','gray'))+
  theme(axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.key.width = unit(1.5,"cm"))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))




##################################################################################
####### OE prediction of different transmission routes on and off ##########
##################################################################################

# When all routes are present (calculated above)

allroutes<-prev_final


## 1. Only Environmental and Adult transfer transmissions present (E + A)

p_h = 0.614
p_v = 0    
lambda = 50 

params_EA = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_EA = dede(func=Monarch_dis,y=ystart,times=times, parms=params_EA)

SA_final_EA = tail(out_EA,n=1)[6] #end of season number of susceptible adults
IA_final_EA = tail(out_EA,n=1)[7] #end of season number of infected adults
CA_final_EA = tail(out_EA,n=1)[8] #end of season number of contaminated adults
NA_final_EA = SA_final_EA + CA_final_EA + IA_final_EA #end of season adult popupaltion size 
prev_final_EA=IA_final_EA/NA_final_EA  # end of season prevalence

EAroutes<-prev_final_EA

## 2. Only Vertical and Adult transfer transmissions present (V + A)

p_h = 0.614
p_v = 0.9
lambda = 0

params_VA = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_VA = dede(func=Monarch_dis,y=ystart,times=times, parms=params_VA)

SA_final_VA = tail(out_VA,n=1)[6] #end of season number of susceptible adults
IA_final_VA = tail(out_VA,n=1)[7] #end of season number of infected adults
CA_final_VA = tail(out_VA,n=1)[8] #end of season number of contaminated adults
NA_final_VA = SA_final_VA + CA_final_VA + IA_final_VA #end of season adult popupaltion size 
prev_final_VA=IA_final_VA/NA_final_VA  # end of season prevalence

VAroutes<-prev_final_VA


## 3. Only Vertical and Environmental transmissions present (V + E)

p_h = 0
p_v = 0.9             
lambda = 50  

params_VE = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_VE = dede(func=Monarch_dis,y=ystart,times=times, parms=params_VE)

SA_final_VE = tail(out_VE,n=1)[6] #end of season number of susceptible adults
IA_final_VE = tail(out_VE,n=1)[7] #end of season number of infected adults
CA_final_VE = tail(out_VE,n=1)[8] #end of season number of contaminated adults
NA_final_VE = SA_final_VE+CA_final_VE+IA_final_VE #end of season adult popupaltion size 
prev_final_VE=IA_final_VE/NA_final_VE  # end of season prevalence

VEroutes<-prev_final_VE


## 4. Only Environmental transmission present (E)

p_h = 0
p_v = 0             
lambda = 50 
params_E = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_E = dede(func=Monarch_dis,y=ystart,times=times, parms=params_E)

SA_final_E = tail(out_E,n=1)[6] #end of season number of susceptible adults
IA_final_E = tail(out_E,n=1)[7] #end of season number of infected adults
CA_final_E = tail(out_E,n=1)[8] #end of season number of contaminated adults
NA_final_E = SA_final_E+CA_final_E+IA_final_E #end of season adult popupaltion size 
prev_final_E=IA_final_E/NA_final_E  # end of season prevalence

Eroutes<-prev_final_E

## 5. Only Adult transfer transmission present (A)

p_h = 0.614
p_v = 0             
lambda = 0  
params_A = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_A = dede(func=Monarch_dis,y=ystart,times=times, parms=params_A)

SA_final_A = tail(out_A,n=1)[6] #end of season number of susceptible adults
IA_final_A = tail(out_A,n=1)[7] #end of season number of infected adults
CA_final_A = tail(out_A,n=1)[8] #end of season number of contaminated adults
NA_final_A = SA_final_A+CA_final_A+IA_final_A #end of season adult popupaltion size 
prev_final_A=IA_final_A/NA_final_A  # end of season prevalence

Aroutes<-prev_final_A

## 6. Only Vertical transmission present (V)

p_h = 0
p_v = 0.9             
lambda = 0  

params_V = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

out_V = dede(func=Monarch_dis,y=ystart,times=times, parms=params_V)

SA_final_V = tail(out_V,n=1)[6] #end of season number of susceptible adults
IA_final_V = tail(out_V,n=1)[7] #end of season number of infected adults
CA_final_V = tail(out_V,n=1)[8] #end of season number of contaminated adults
NA_final_V = SA_final_V + CA_final_V + IA_final_V #end of season adult popupaltion size 
prev_final_V=IA_final_V/NA_final_V  # end of season prevalence

Vroutes<-prev_final_V

#Make a dataframe with route combinations
require(reshape2)
routes<-data.frame(allroutes,EAroutes,VAroutes,VEroutes,Eroutes, Aroutes,Vroutes)

#Format for ggplot
routes<-melt(routes)
routes <- routes[order(-routes$value),] #ordered highest to lowest prevalence


#Figure 3b
ggplot()+
  geom_bar(data=routes,aes(x=variable, y=value),na.rm=TRUE, stat="identity")+
  theme_classic()+
  labs(fill="", x="Transmission route present", y = "Proportion infected")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 0.8))+
  theme(axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=19, hjust = 0.5),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.key.width = unit(1.5,"cm"))+
  scale_x_discrete(labels=c("All","E+A", "V+A", "V+E","E","A","V" ))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))



######################################################################################################
## Virulence effects on adult abundance and prevalence with different combos of transmission routes ##
######################################################################################################
## Save predictions with all transmission routes ##
#########################################################

testVad_allroutes = NULL

p_h = 0.614
p_v = 0.9
lambda = 50

for (z in seq(1,24,1)){ # vary virulence via infected adult mortality
   mu_I = 1/z
   params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

   ystart = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=18,I_A=2,C_A=0,M=7000,W=0)
   
  out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)

  time<-out[,1]
  SA_final = tail(out,n=1)[6]
  CA_final = tail(out,n=1)[8]
  IA_final = tail(out,n=1)[7]

# calculate total adult pop size at end of season
  NA_final = SA_final+CA_final+IA_final
  
  # calculate prevalence
  prev_final = IA_final/NA_final
  prev_final
  
  # calculate virulence
  virulence = 1 - ((b_i/b_s) * (mu_S/mu_I) * theta)

  testVad_allroutes=rbind(testVad_allroutes,c(z, virulence, NA_final, prev_final))

      }

#save output
tnAll<-as.data.frame(testVad_allroutes)
colnames(tnAll) <- c("mu_I", "virulence","Abundance", "Prevalence")


#########################################################
## Save predictions with no adult transfer transmission##
#########################################################


testVad_noph = NULL


for (k in seq(1,24,1)){
  mu_I = 1/k
  p_h=0
  params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)
  
  ystart = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=18,I_A=2,C_A=0,M=7000,W=0)
  
  out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)

  SA_final = tail(out,n=1)[6]
  CA_final = tail(out,n=1)[8]
  IA_final = tail(out,n=1)[7]

  NA_final = SA_final+CA_final+IA_final
  
  # calculate prevalence
  prev_final = IA_final/NA_final
  
  # calculate virulence
  virulence = 1 - ((b_i/b_s) * (mu_S/mu_I) * theta)  
  
  testVad_noph=rbind(testVad_noph,c(k, virulence, NA_final, prev_final))
  
}

tnnoph<-as.data.frame(testVad_noph)
colnames(tnnoph) <- c("mu_I", "virulence","Abundance", "Prevalence")


#########################################################
## Save predictions with no environmental transmission###
#########################################################

testVad_nope = NULL

for (l in seq(1,24,1)){
  mu_I = 1/l
  p_h = 0.614
  lambda=0
  # define the vector of parameters
  params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)
  
  ystart = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=18,I_A=2,C_A=0,M=7000,W=0)
  
  out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
  

  SA_final = tail(out,n=1)[6]
  CA_final = tail(out,n=1)[8]
  IA_final = tail(out,n=1)[7]

  NA_final = SA_final+CA_final+IA_final
  
  # calculate prevalence
  prev_final = IA_final/NA_final
  
  # calculate virulence
  virulence = 1 - ((b_i/b_s) * (mu_S/mu_I) * theta)  
  
  testVad_nope=rbind(testVad_nope,c(l, virulence, NA_final, prev_final))
  
}

tnnope<-as.data.frame(testVad_nope)
colnames(tnnope) <- c("mu_I", "virulence","Abundance", "Prevalence")

# Assign names to vectors with pedictions
tnnope$route <- rep("V + A", 24) # no environmental transmission
tnnoph$route <- rep("V + E", 24) # no adult transfer  
tnAll$route <- rep("All routes", 24) # all routes present

#Combine the three vectors with predictions
tn<-rbind(tnnope,tnnoph)
tn<-rbind(tn,tnAll)



# Figure 4
p1<- ggplot(tn, aes(x=virulence, y=Abundance/NA_final_f))+ # NA_final_f is disease-free population size 
  geom_line(aes(color = route),size = 2,na.rm=TRUE)+
  geom_vline(xintercept = 0.58, color = 'gray',linetype="dashed")+
labs(fill="", x="Virulence", y = "Relative host abundance")+
  theme_classic()+
  scale_color_manual(values=c('black','cyan4', 'orange'))+
  theme(axis.text=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=16, hjust = 1),
        axis.title.y = element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+ ##top, right, bottom, left
  scale_y_continuous(limits = c(0.3, 1))

p2<- ggplot(tn, aes(x=virulence, y=Prevalence))+
  geom_line(aes(color=route),size = 2,na.rm=TRUE)+
  geom_vline(xintercept = 0.58, color = 'gray',linetype="dashed")+
labs(fill="", x="Virulence", y = "Proportion infected")+
  theme_classic()+
  scale_color_manual(values=c('black','cyan4', 'orange'))+
  theme(axis.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=16, hjust = 1),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=14),
        legend.key.width = unit(.5,"cm"),
        legend.title=element_blank(),
        legend.position = c(.7,.9),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  scale_y_continuous(limits = c(0, 1))
grid.arrange(p1, p2, nrow = 1)



#################################################

###### FIGURES FOR SUPPLEMENTAL MATERIALS ######

#################################################

###clear all workspace again
rm(list=ls()) 
graphics.off()

###load necessary packages
require(deSolve)

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
  M = y[8] # milkweed
  W = y[9] # contaminated milkweed
  
  
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
larvdens = 1.3          # Larval density
rho = 0.76              # Probability of pupa surviving to adult
mu_S = 1/24             # Mortality of uninfected adult
mu_I = 1/20             # Mortality rate of infected adult
mu_d = -log(0.06)/9     # Density-independent per capita larval mortality rate 
mu_1 = (200/larvdens)*((b_s*p_se*g*rho/mu_S) - mu_d - g)   
# Density-dependent per capita larval mortality rate 
tau2=7                  # Time in pupa stage
theta = 0.72            # Probability of infected adult eclosing and mating
p_v = 0.9               # Probability of vertical transmission
p_h =  0.614            # Probability of contaminated adult infecting their offspring
delta = 1/2             # Daily adult mating probability  
mu_C= 1/14              # Rate of spore loss for contaminated adults
r = 0.032               # Milkweed growth rate
K= 40000                # End of season number of milkweed leaves (assuming 200 leaves per plant)
c = 35/9                # Larval consumption rate of milkweed
lambda = 50             # Milkweed leaf visitation rate by infected adults,  (leaves/day)
mu_W = 0.0125           # Spore decay rate on milkweed 


# Define the range of time over which we want to solve the equations
times = seq(0,150,by=1) 

# Initial conditions
ystart = c(S_L=0,I_V=0,I_E=0,I_H=0,S_A=18,I_A=2,C_A=0,M=7000,W=0)

# Combine paramameters into single vector
params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)


#Figures S1

####################################################################################################################
###Vary Duration of adult contamination (1/mu_C) and Probability of contaminated adult infecting offspring (ph) ####
####################################################################################################################
# let z be mu_C
z=seq(1,15,length.out = 25)

# let q be p_h 
q=seq(0,1,length.out = 25)


# create an empty matrix to store results
total=matrix(rep(0,length(q)*length(z)),nrow=length(q)) #prevalence
adults=matrix(rep(0,length(q)*length(z)),nrow=length(q)) #abundance

for(i in 1:length(q)){ # vary p_h
  p_h =q[i]
  for(j in 1:length(z)){# vary mu_C
    mu_C = 1/z[j]
    
    # define parameters
    params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

    out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
    
    SA_f = tail(out,n=1)[6]
    CA_f = tail(out,n=1)[8]
    IA_f = tail(out,n=1)[7]
    
    # calculate total adult population size 
    NA_f = SA_f+CA_f+IA_f

    # write prevalence to i,j entry of matrix
    total[i,j] = IA_f/NA_f #prevalence
    adults[i,j] = NA_f #abundance
  } # end j loop
} # end i loop



#FigureS1a
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right

NA_final_f = 164.6708 

filled.contour(q,z,adults/NA_final_f,
               nlevels=15,
               col=rev(heat.colors(20)),
               plot.title=title(main = "Relative adult abundance",
                                xlab=expression(paste("Probability of contaminated ")),
                                ylab=expression(paste("Duration of adult contamination  ", (T[c]))),
                                cex.main=1.2,cex.lab=1.2))
mtext(expression(paste("    adult infecting offspring ",(p[h]))),side = 1, line = 4, adj = 0, cex = 1.2)



#FigureS1b
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(q,z,total,
               nlevels=15,
               col=rev(heat.colors(15)),
               plot.title=title(main = "Proportion infected adults",
                                xlab=expression(paste("Probability of contaminated")),
                                ylab=expression(paste("Duration of adult contamination   ", (T[c]))),
                                cex.main=1.2,cex.lab=1.2))
mtext(expression(paste("    adult infecting offspring ",(p[h]))),side = 1, line = 4, adj = 0, cex = 1.2)



########################################################################################
###Milkweed visitation rate and Duration of milkweed contamination (1/mu_W) ####
########################################################################################
# let w be lambda
w=seq(10,100,length.out = 25) 

# let x be mu_w
x=seq(1,80,length.out = 25)

# create an empty matrix to store adult prevalence (ad_prev)
tot3=matrix(rep(0,length(w)*length(x)),nrow=length(w))
ad3=matrix(rep(0,length(w)*length(x)),nrow=length(w))

for(i in 1:length(w)){ # vary 
  lambda =w[i]
  for(j in 1:length(x)){# vary 
    mu_W = 1/x[j]
    # define parameters
    params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)
    # solve def eq
    
    out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
    
    SA_final = tail(out,n=1)[6]
    CA_final = tail(out,n=1)[8]
    IA_final = tail(out,n=1)[7]
    NA_final = SA_final+CA_final+IA_final
    
    # calculate prevalence
    prev_final=IA_final/NA_final  
    
    # write prevalence to i,j entry of matrix
    tot3[i,j] = IA_final/NA_final
    ad3[i,j] = NA_final
  } # end j loop
} # end i loop

#FigureS1c
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3.1,3))# bottom, left, top, right
filled.contour(w,x,ad3/NA_final_f,
               nlevels=15,
               col=rev(heat.colors(20)),
               plot.title=title(main = "Relative adult abundance",
                                xlab=expression(paste("Milkweed visitation rate ",(lambda))),
                                ylab=expression(paste("Duration of milkweed contamination ", (T[w]))),
                                cex.main=1.2,cex.lab=1.2))



#FigureS1d
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3.1,3))# bottom, left, top, right
filled.contour(w,x,tot3,
               nlevels=10,
               col=rev(heat.colors(15)),
               plot.title=title(main = "Proportion infected adults",
                                xlab=expression(paste("Milkweed visitation rate ",(lambda))),
                                ylab=expression(paste("Duration of milkweed contamination ", (T[w]))),
                                cex.main=1.2,cex.lab=1.2))

#################################################################################################
###  Probability of contaminated adult infecting offspring (ph) and Milkweed visitation rate ####
#################################################################################################
# let d be ph
d=seq(0,1,length.out = 25)

# let v be lambda
v=seq(10,100,length.out = 25)

adul2=matrix(rep(0,length(d)*length(v)),nrow=length(d))
tot2=matrix(rep(0,length(d)*length(v)),nrow=length(d))


for(i in 1:length(d)){ # vary ph
  p_h = d[i]
  for(j in 1:length(v)){# vary lambda
    lambda = v[j]
    
    # define parameters
    
    params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)  
    
    out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
    
    SA_final = tail(out,n=1)[6]
    CA_final = tail(out,n=1)[8]
    IA_final = tail(out,n=1)[7]
    
    # calculate total adult pop size at t=150
    NA_final = SA_final+CA_final+IA_final
    
    # calculate prevalence
    prev_final=IA_final/NA_final  
    
    # write prevalence to i,j entry of matrix
    tot2[i,j] = IA_final/NA_final
    
    adul2[i,j] = NA_final
  } # end j loop
} # end i loop


#Figure_S1e
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(d,v,tot2,
               nlevels=15,
               col=rev(heat.colors(20)),
               plot.title=title(
                 ylab=expression(paste("Milkweed visitation rate ", lambda)),
                 main = "Proportion infected adults",
                 xlab = expression(paste("Probability of contaminated")),
                 cex.main=1.2,cex.lab=1.2))
mtext(expression(paste("    adult infecting offspring ",(p[h]))),side = 1, line = 4, adj = 0, cex = 1.2)


#Figure_S1f
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(d,v,adul2/NA_final_f,
               nlevels=10,
               col=rev(heat.colors(15)),
               plot.title=title(ylab=expression(paste("Milkweed visitation rate ", lambda)),
                                main = "Relative adult abundance",
                                xlab = expression(paste("Probability of contaminated adult")),
                                cex.main=1.2,cex.lab=1.2))
mtext(expression(paste("    adult infecting offspring ",(p[h]))),side = 1, line = 4, adj = 0, cex = 1.2)

####################################################################################################################

#Figures S2

###################################################################################
###  Probability of contaminated adult infecting offspring (ph) and Virulence ####
###################################################################################

# let h be ph
h=seq(0,1,length.out = 25)

# let p be mu_I (virulence component)
p=seq(from = 1/24,to = 1,by= 0.01) 

ad4=matrix(rep(0,length(p)*length(h)),nrow=length(p))
to4=matrix(rep(0,length(p)*length(h)),nrow=length(p))
virul=NULL

for(i in 1:length(h)){ # vary ph 
  p_h = h[i]
  for(j in 1:length(p)){# vary virulence 
    mu_I = p[j]
    
    params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

    out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
    
    SA_final = tail(out,n=1)[6]
    CA_final = tail(out,n=1)[8]
    IA_final = tail(out,n=1)[7]
    
    NA_final = SA_final+CA_final+IA_final
    
    # calculate prevalence
    prev_final=IA_final/NA_final  
    
    # write results to i,j entry of matrix
    to4[j, i] = IA_final/NA_final
    ad4[j,i] = NA_final
    virul[j] <- 1 - ((b_i/b_s) * (mu_S/mu_I) * theta)
  } # end j loop
} # end i loop

#Figure S2a
par(oma=c(0,1,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(x= virul, y = h, ad4/NA_final_f,
               nlevels=20,
               col=rev(heat.colors(20)),
               plot.title=title(main =  "Relative adult abundance", 
                                ylab=expression(paste("adult infecting offspring ",(p[h]))),
                                xlab="Virulence",
                                cex.main=1.2,cex.lab=1.2),
               key.axes = axis(4, seq(0.2, 1, by = 0.10),asp=1))
mtext(expression(paste("Probability of contaminated")),side = 2, line = 4, adj = 0, cex = 1.2)


#Figure S2b
par(oma=c(0,1,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(y = h, x = virul, to4,
               nlevels=20,
               col=rev(heat.colors(20)),
               plot.title=title(main = "Proportion infected adults", 
                                ylab=expression(paste("adult infecting offspring ",(p[h]))),
                                xlab="Virulence",
                                cex.main=1.2,cex.lab=1.2))
mtext(expression(paste("Probability of contaminated")),side = 2, line = 4, adj = 0, cex = 1.2)



################################################
###  Milkweed visitation rate and Virulence ####
################################################


# let n be lambda
n=seq(10,100,length.out = 25) 

#let o be mu_I (virulence component)
o=seq(from = 1/24,to = 1,length.out = 25) 

# create an empty matrix to store results

adults5=matrix(rep(0,length(o)*length(n)),nrow=length(o))
total5=matrix(rep(0,length(o)*length(n)),nrow=length(o))
virul2=NULL

for(i in 1:length(n)){ # vary lamd 
  lambda = n[i]
  for(j in 1:length(o)){# vary virulence 
    mu_I = o[j]
    
    # define parameters
    params = c(r, K, b_s,b_i,mu_d,mu_1,g,c,theta,lambda,p_v,mu_I,mu_S,mu_W,mu_C,delta,p_h,tau1,tau2)

    out = dede(func=Monarch_dis,y=ystart,times=times, parms=params)
    
    SA_final = tail(out,n=1)[6]
    CA_final = tail(out,n=1)[8]
    IA_final = tail(out,n=1)[7]
    
    NA_final = SA_final+CA_final+IA_final
    
    # calculate prevalence
    prev_final=IA_final/NA_final  
    
    
    # write results to i,j entry of matrix
    total5[j, i] = IA_final/NA_final
    adults5[j,i] = NA_final
    virul2[j] <- 1 - ((b_i/b_s) * (mu_S/mu_I) * theta)
  } # end j loop
} # end i loop


#Figure S2c
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(x= virul2, y = n, adults5/NA_final_f,
               nlevels=20,
               col=rev(heat.colors(20)),
               plot.title=title(main =  "Relative adult abundance", 
                                ylab=expression(paste("Milkweed visitation rate ", lambda)),
                                xlab="Virulence ",
                                cex.main=1.2,cex.lab=1.2))
#Figure S2d
par(oma=c(0,0.25,0,0), mar=c(5,4.25,3,3))# bottom, left, top, right
filled.contour(y = n, x = virul2, total5,
               nlevels=20,
               col=rev(heat.colors(20)),
               plot.title=title(main = "Proportion infected adults", ylab=expression(paste("Milkweed visitation rate ", lambda)),
                                xlab="Virulence",
                                cex.main=1.2,cex.lab=1.2))




