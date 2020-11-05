#Code for all analyses and the one figure in :
#"The effects of exposure to predators on personality and plasticity"
#Amy Bucklaew & Ned A. Dochtermann

#Requires data=BucklaewData2018, available as "BucklaewData2018.csv"
#
library(MCMCglmm); library(lme4); library(rptR);library(lmerTest)

#### Misc Functions ####
#Calculating repeatabilities with fixed effects variances
#and distribution specific variances in denominator.
#The functions below calculate these variances.
fixed.var<-function(mod=model){
  vmVarF<-numeric(1000)
  for(i in 1:1000){
    Var<-var(as.vector(mod$Sol[i,] %*% t(mod$X)))
    vmVarF[i]<-Var}
  as.mcmc(vmVarF)
}


Pois.win.var<-function(mod=model,with.var=2){
  dist.var<-mod$Sol[,1]
  dist.var<-log((1/exp(dist.var))+1)
  omega<-mod$VCV[,with.var]
  tot.win.var<-as.mcmc(omega+dist.var)
}

Bin.win.var <- function(mod=model,with.var=2){
  dist.var<-(pi^2)/3
  omega<-mod$VCV[,with.var]
  tot.win.var<-as.mcmc(omega+dist.var)
}

#### Scale fixed effects & reorder factors so "Pre" is the index factor ####
BucklaewData2018$TempSCALED <- scale(BucklaewData2018$Temp)
BucklaewData2018$MassSCALED <- scale(BucklaewData2018$Mass)
BucklaewData2018$Exp <- as.factor(BucklaewData2018$Exp)
BucklaewData2018$Exp <- factor(BucklaewData2018$Exp,levels(BucklaewData2018$Exp)[c(2,1)])

#### MCMC evaluation of difference in repeatability by behav type ####
#MCMC Chain & Prior Specification
NITT=13000;THIN=10;BURN=3000
multi=200
a=10000

prior.uni<-list(G=list(G1=list(V=1,nu=1, alpha.mu=0, alpha.V=a)),
                R=list(V=1,nu=0.002))

prior.bi<-list(G=list(G1=list(V=diag(2),nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)),
               R=list(V=diag(2),nu=1.002))

prior.3<-list(G=list(G1=list(V=diag(2),nu=1.002, alpha.mu=c(0,0), alpha.V=diag(2)*a)),
              R=list(V=1,nu=0.002))

prior.4<-list(G=list(G1=list(V=1,nu=0.002, alpha.mu=0, alpha.V=a)),
              R=list(V=diag(2),nu=1.002))

#Distance traveled (Activity)
#Model fitting
Dist1.mcmc<-MCMCglmm(Dist~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                     random=~ID,rcov=~units,
                     nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                     family="gaussian",prior=prior.uni,
                     verbose=F,data=BucklaewData2018)

Dist2.mcmc<-MCMCglmm(Dist~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                     random=~idh(Exp):ID,rcov=~units,
                     nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                     family="gaussian",prior=prior.3,
                     verbose=F,data=BucklaewData2018)

Dist3.mcmc<-MCMCglmm(Dist~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                     random=~ID,rcov=~idh(Exp):units,
                     nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                     family="gaussian",prior=prior.4,
                     verbose=F,data=BucklaewData2018)

Dist4.mcmc<-MCMCglmm(Dist~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                     random=~idh(Exp):ID,rcov=~idh(Exp):units,
                     nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                     family="gaussian",prior=prior.bi,
                     verbose=F,data=BucklaewData2018)

#Model Comparison 
round(Dist1.mcmc$DIC,2)
round(Dist2.mcmc$DIC,2)
round(Dist3.mcmc$DIC,2)
round(Dist4.mcmc$DIC,2)

round(Dist1.mcmc$DIC,2)-round(Dist3.mcmc$DIC,2)
round(Dist2.mcmc$DIC,2)-round(Dist3.mcmc$DIC,2)
round(Dist3.mcmc$DIC,2)-round(Dist3.mcmc$DIC,2)
round(Dist4.mcmc$DIC,2)-round(Dist3.mcmc$DIC,2)

#Dist Variance Components & Repeatability: Estimate from Dist4.mcmc
Vi.Dist<-Dist4.mcmc$VCV[,c(1,2)]
Vf.Dist<-fixed.var(mod=Dist4.mcmc)
Vw.Dist<-Dist4.mcmc$VCV[,c(3,4)]

posterior.mode(Vi.Dist)
posterior.mode(Vf.Dist)
posterior.mode(Vw.Dist)
posterior.mode(Vi.Dist/(Vi.Dist+Vf.Dist+Vw.Dist))

HPDinterval(Vi.Dist)
HPDinterval(Vf.Dist)
HPDinterval(Vw.Dist)
HPDinterval(Vi.Dist/(Vi.Dist+Vf.Dist+Vw.Dist))

#Unique Zones
#Model Fitting
UZ1.mcmc<-MCMCglmm(UZ~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                   random=~ID,rcov=~units,
                   nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                   family="poisson",prior=prior.uni,
                   verbose=F,data=BucklaewData2018)

UZ2.mcmc<-MCMCglmm(UZ~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                   random=~idh(Exp):ID,rcov=~units,
                   nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                   family="poisson",prior=prior.3,
                   verbose=F,data=BucklaewData2018)

UZ3.mcmc<-MCMCglmm(UZ~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                   random=~ID,rcov=~idh(Exp):units,
                   nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                   family="poisson",prior=prior.4,
                   verbose=F,data=BucklaewData2018)

UZ4.mcmc<-MCMCglmm(UZ~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                   random=~idh(Exp):ID,rcov=~idh(Exp):units,
                   nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                   family="poisson",prior=prior.bi,
                   verbose=F,data=BucklaewData2018)

round(UZ1.mcmc$DIC,2)
round(UZ2.mcmc$DIC,2)
round(UZ3.mcmc$DIC,2)
round(UZ4.mcmc$DIC,2)

round(UZ1.mcmc$DIC,2)-round(UZ2.mcmc$DIC,2)
round(UZ2.mcmc$DIC,2)-round(UZ2.mcmc$DIC,2)
round(UZ3.mcmc$DIC,2)-round(UZ2.mcmc$DIC,2)
round(UZ4.mcmc$DIC,2)-round(UZ2.mcmc$DIC,2)

#UZ Variance Components & Repeatability: Estimate from UZ4.mcmc
Vi.UZ<-UZ4.mcmc$VCV[,c(1,2)]
Vf.UZ<-fixed.var(mod=UZ4.mcmc)
Vw.UZ<-Pois.win.var(mod=UZ4.mcmc,with.var=c(3,4))

posterior.mode(Vi.UZ)
posterior.mode(Vf.UZ)
posterior.mode(Vw.UZ)
posterior.mode(Vi.UZ/(Vi.UZ+Vf.UZ+Vw.UZ))

HPDinterval(Vi.UZ)
HPDinterval(Vf.UZ)
HPDinterval(Vw.UZ)
HPDinterval(Vi.UZ/(Vi.UZ+Vf.UZ+Vw.UZ))


#Latency
#Model Fitting
Lat1.mcmc<-MCMCglmm(LatBi~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                    random=~ID,rcov=~units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    family="categorical",prior=prior.uni,
                    verbose=F,data=BucklaewData2018)

Lat2.mcmc<-MCMCglmm(LatBi~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                    random=~idh(Exp):ID,rcov=~units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    family="categorical",prior=prior.3,
                    verbose=F,data=BucklaewData2018)

Lat3.mcmc<-MCMCglmm(LatBi~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                    random=~ID,rcov=~idh(Exp):units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    family="categorical",prior=prior.4,
                    verbose=F,data=BucklaewData2018)

Lat4.mcmc<-MCMCglmm(LatBi~Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED,
                    random=~idh(Exp):ID,rcov=~idh(Exp):units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    family="categorical",prior=prior.bi,
                    verbose=F,data=BucklaewData2018)

#Model Comparison
round(Lat1.mcmc$DIC,2)
round(Lat2.mcmc$DIC,2)
round(Lat3.mcmc$DIC,2)
round(Lat4.mcmc$DIC,2)

round(Lat1.mcmc$DIC,2)-round(Lat2.mcmc$DIC,2)
round(Lat2.mcmc$DIC,2)-round(Lat2.mcmc$DIC,2)
round(Lat3.mcmc$DIC,2)-round(Lat2.mcmc$DIC,2)
round(Lat4.mcmc$DIC,2)-round(Lat2.mcmc$DIC,2)

#Lat Variance Components & Repeatability: Estimate from Lat4.mcmc
Vi.Lat<-Lat4.mcmc$VCV[,c(1,2)]
Vf.Lat<-fixed.var(mod=Lat4.mcmc)
Vw.Lat<-Bin.win.var(mod=Lat4.mcmc,with.var=c(3,4))

posterior.mode(Vi.Lat)
posterior.mode(Vf.Lat)
posterior.mode(Vw.Lat)
posterior.mode(Vi.Lat/(Vi.Lat+Vf.Lat+Vw.Lat))

HPDinterval(Vi.Lat)
HPDinterval(Vf.Lat)
HPDinterval(Vw.Lat)
HPDinterval(Vi.Lat/(Vi.Lat+Vf.Lat+Vw.Lat))

#### Plotting ####

#create tables for plotting
Dist.ef <- matrix(c(posterior.mode(Dist3.mcmc$Sol[,1]),posterior.mode(Dist3.mcmc$Sol[,1]+Dist3.mcmc$Sol[,2]),
                    median(Dist3.mcmc$Sol[,1]),median(Dist3.mcmc$Sol[,1]+Dist3.mcmc$Sol[,2]),
                    HPDinterval(Dist3.mcmc$Sol[,1])[1],HPDinterval(Dist3.mcmc$Sol[,1]+Dist3.mcmc$Sol[,2])[1],            
                    HPDinterval(Dist3.mcmc$Sol[,1])[2],HPDinterval(Dist3.mcmc$Sol[,1]+Dist3.mcmc$Sol[,2])[2]),2,4) 

UZ.ef <- matrix(c(posterior.mode(UZ2.mcmc$Sol[,1]),posterior.mode(UZ2.mcmc$Sol[,1]+UZ2.mcmc$Sol[,2]),
                    median(UZ2.mcmc$Sol[,1]),median(UZ2.mcmc$Sol[,1]+UZ2.mcmc$Sol[,2]),
                    HPDinterval(UZ2.mcmc$Sol[,1])[1],HPDinterval(UZ2.mcmc$Sol[,1]+UZ2.mcmc$Sol[,2])[1],            
                    HPDinterval(UZ2.mcmc$Sol[,1])[2],HPDinterval(UZ2.mcmc$Sol[,1]+UZ2.mcmc$Sol[,2])[2]),2,4) 


#because MCMCglmm fits fixed effects in a different way for categorical, 
#coefficients aren't directly convertable using plogis
#instead, plot effects and 95%CI's with effects package from a glmer fitted model
#(direction and significance of effects are the same, regardless of fitting method)
m.lat <- glmer(LatBi~-1+Exp+TempSCALED+MassSCALED+Exp:TempSCALED+Exp:MassSCALED+(Exp|ID),family=binomial,data=BucklaewData2018)
lat.95 <- confint(m.lat,parm="beta_")
lat.ef <- plogis(cbind(fixef(m.lat)[1:2],lat.95[1:2,]))#convert to probability

#Figure 1
#tiff('F1.tiff',width=8,height=10,units='in',res=600)
par(mfrow=c(3,1),pty='s',mar=c(1,1,1,1)+.01,oma=c(3,0,0,0))
plot(NA, xaxt="n",xlab="",ylab="Distance (cm)",
     xlim=c(0.75,2.25), ylim=c(floor(min(Dist.ef[,4])),ceiling(max(Dist.ef[,5]))),
     cex.axis=1.25,cex.lab=1.5)
arrows(1:2,Dist.ef[,4],1:2,Dist.ef[,5],length=.1,angle=90,code=3)
lines(c(1,2),Dist.ef[,2],lty=2)
points(c(1,2),Dist.ef[,2],cex=3,pch=21,bg='white')
axis(1, at=c(1,2),labels=F)
mtext("A",side=2,line=4.5,padj=-7,cex=1.5,las=2)

plot(NA, xaxt="n",xlab="",ylab="Unique Zones Visited",
     xlim=c(0.75,2.25),ylim=c(floor(min(UZ.ef[,4])),ceiling(max(UZ.ef[,5]))),
     cex.axis=1.25,cex.lab=1.5)
arrows(1:2,UZ.ef[,4],1:2,UZ.ef[,5],length=.1,angle=90,code=3)
lines(c(1,2),UZ.ef[,2],lty=2)
points(c(1,2),UZ.ef[,2],cex=3,pch=21,bg='white')
axis(1, at=c(1,2),labels=F)
mtext("B",side=2,line=4.5,padj=-7,cex=1.5,las=2)

plot(NA, xaxt="n",xlab="",ylab="Emergence (proportion)",
     xlim=c(0.75,2.25),ylim=c(0,1),
     cex.axis=1.25,cex.lab=1.5)
arrows(1:2,lat.ef[,2],1:2,lat.ef[,3],length=.1,angle=90,code=3)
lines(c(1,2),lat.ef[,1],lty=2)
points(c(1,2),lat.ef[,1],cex=3,pch=21,bg='white')
axis(1, at=c(1,2),labels=c("Pre-Exposure","Post-Exposure"),cex.axis=1.5)
mtext("C",side=2,line=4.5,padj=-7,cex=1.5,las=2)
#dev.off()

#### Supplemental post-hoc analysis (see Response to Reviewers) ####
# Calculate delta-V per:
# Royauté, R., and N. A. Dochtermann. 2020. 
# Comparing ecological and evolutionary variability within datasets.  
# ecoevorxiv.org doi:10.32942/osf.io/tn7u5
# https://ecoevorxiv.org/tn7u5/

#Activity
Dist_Vi_diff <- Dist4.mcmc$VCV[,1]-Dist4.mcmc$VCV[,2]
Dist_Vw_diff <- Dist4.mcmc$VCV[,3]-Dist4.mcmc$VCV[,4]
posterior.mode(Dist_Vi_diff); posterior.mode(Dist_Vw_diff)
sum(Dist_Vi_diff<0); sum(Dist_Vw_diff<0)
#Model comparison demonstrated support for Vi & Vw differing between
#pre- and post-exposure. This was supported in 988 of 1000 posterior estimates for Vw (very strong support)
#and 884/1000 for Vi


#Unique Zones
UZ_Vi_diff <- UZ4.mcmc$VCV[,1]-UZ4.mcmc$VCV[,2]
UZ_Vw_diff <- UZ4.mcmc$VCV[,3]-UZ4.mcmc$VCV[,4]
posterior.mode(UZ_Vi_diff); posterior.mode(UZ_Vw_diff)
sum(UZ_Vi_diff<0); sum(UZ_Vw_diff<0)
#Model comparison demonstrated some limited support for Vw differeing between
#pre- and post-exposure. This difference was supported in 903 of 1000 posterior estimates (moderate support)

#Latency
Lat_Vi_diff <- Lat4.mcmc$VCV[,1]-Lat4.mcmc$VCV[,2]
Lat_Vw_diff <- Lat4.mcmc$VCV[,3]-Lat4.mcmc$VCV[,4]
posterior.mode(Lat_Vi_diff); posterior.mode(Lat_Vw_diff)
sum(Lat_Vi_diff<0); sum(Lat_Vw_diff>0)
#Model comparison demonstrated weak support for Vi differeing between
#pre- and post-exposure. This difference was supported in 905 of 1000 posterior estimates (moderate support)
