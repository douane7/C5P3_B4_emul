# C5P3_B4_emul6: based on raw CMIP5/PMIP3+Biome4 simulations 
# (from files provides by CMIP5_biomes3.R and PMIP3_Biomes3.R and PMIP3_1k_biomes3.R)
# aggregate the biome output to be used to calaculate the emulator (geospatial regression GWR)
# - based random taking from the gridded data by biome
# input variables: present clim, CONCENTRATIONS, orbital parameters, pop
# output variables: bioclimate, NPP 
# 04/08/2020
rm(list = ls())
require(gdata)
require(fields)
library(raster)
library(maptools)
library( GWmodel )
library(akima)

source("cartosurf.r")
source('assignbiome2.r')
RDfile="C5P3_B4_GWRall.RData"    # used only to receive all the calculated data
RDfil2="C5P3_B4_GWRmod5.RData"
RDtxt="C5P3_B4_GWRall5.txt"
RDpdf="C5P3_B4_GWRall5.pdf"
pdf(RDpdf)
cat(RDtxt,'\n',file=RDtxt)
capture.output(date(),file=RDtxt,append=T)
Driv2=read.table('driv-proj6.txt',header=T,row.names=1)
np=nrow(Driv2)
win= c(-20, 60, -5, 75)

# reading forcing data (GHG concentration) Present future
x = read.xls ("R26_bulk.xls", sheet = 1, header = TRUE)
F26=x[c(1,2,3),c(5:16)]
rownames(F26)=c("CO2","CH4","N2O")
x = read.xls ("R45_bulk.xls", sheet = 1, header = TRUE)
F45=x[c(1,2,3),c(5:16)]
rownames(F45)=rownames(F26)
x = read.xls ("R85_bulk.xls", sheet = 1, header = TRUE)
F85=x[c(1,2,3),c(5:16)]
rownames(F85)=rownames(F26)

# global population  2000, 2005,2010, 2020, ... for each RCP
#x 1000 hab (Van Vuuren et al 2011)  green/red/blue
pop=matrix(NA,3,12)
pop[,1]=c(6144834,6144834,6144834)  # 2000
pop[,2]=c(6542159,6542159,6542159)  # 2005
pop[,3]=c(6958169,6958169,6958169)  # 2010
pop[,4]=c(7510000,7505000,7530000)  # 2020
pop[,5]=c(8200000,8180000,8800000)  # 2030
pop[,6]=c(8800000,8500000,9400000)  # 2040
pop[,7]=c(9000000,8900000,10400000)  # 2050
pop[,8]=c(9100000,9000000,10700000)  # 2060
pop[,9]=c(9150000,9050000,11500000)  # 2070
pop[,10]=c(9150000,9050000,11900000)  # 2080
pop[,11]=c(9150000,8950000,12050000)  # 2090
pop[,12]=c(9050000,8800000,12300000)  # 2100
colnames(pop)=c(2000,2005,seq(2010,2100,10))
rownames(pop)=c("RCP26","RCP45","RCP85")


# paleo (lgm,midHol,1K)
Fpal=matrix(NA,7,6)     # 
Fpal[,1]=c(185,280,280,280,280,280,280)  #CO2
Fpal[,2]=c(350,650,650,650,650,650,650)  #CH4
Fpal[,3]=c(200,270,270,270,270,270,270)  #N2O 
Fpal[,4]=c(2000,28000,394250,395778,389774,767975,1082361)  # population Hyde 3.1
Fpal[,5]=c(0.3,0.3,0.2,0.4,0.2,0.2,0.5)  # Volcanoes (midHol= moyenne 1960-2000)
Fpal[,6]=c(270,270,400,200,200,600,400)  # solar (midHol= moyenne 1960-2000)
#Fpal[,6]=c(0,0,-0.5,-0.9,0.1,-0.3)  # solar Steinhilber 2012 anomalies W/m2
rownames(Fpal)=c("lgm","midHol","1200","1255","1380","1720","1810")
colnames(Fpal)=c(rownames(F26),'POP','Volc','Solar')


# orbital parameters 21k, 6k, 0k
Ins=matrix(c(0.0196,0.0010,23.4,-213,                #19k
             0.0190,-0.005,24.2,-27,                #7k
             0.0170,0.0170,23.6,85,                #1k
             0.0170,0.0170,23.6,85,                #1k
             0.0170,0.0170,23.6,85,                #1k
             0.0170,0.0170,23.6,85,                #1k
             0.0167,0.0164,23.4,102),4,7)          #0k
rownames(Ins)=c("Ecc","Prec","Obl","Omega")
colnames(Ins)=c("lgm","midHol","V+S-","V-S-","V-S+","V+S+","Present")
iins=c(1,3,4)   # we do not use precession (Prec)

#########################################
# using biome4  (modern) 
PFTNAMES=c('tet','trt','tbe','tst','ctc','bec','bst','tg','trg4','wds','tsg','ch','lf')
APFTNAMES=c('trt','tbe','tet','bot','teg','trg4','wds','tug')
b4names=c('TrEgFo','TrSeDeFo','TrDeFo','TeDeFo','TeCoFo',
          'WaMxFo','CoMxFo','CoCoFo','ClMxFo','BoEgFo',
          'BoDeFo','TrSav','TrXeShl', 'TeXsShl','TeXsWol',
          'TeBlSav','CoOpWol','TeDes','TrGrl','TeGrl',
          'HotDes','TundSt','TundSh','TundDSh','TundPsSh',
          'TundCu', 'BarLand','Ice')
codeb=matrix(NA,28,2)
codeb[,1]=c('darkslategrey','darkolivegreen','yellowgreen','lightgreen','chartreuse3',
            'blue4','dodgerblue','cyan3','thistle3','slateblue3',
            'slateblue1','darkred','darkgoldenrod4','darkorange3','brown2',
            'brown1','thistle','khaki1','burlywood3','burlywood1',
            'cornsilk','magenta','darkorchid4','darkorchid3','darkorchid1',
            'darkorchid','gray86','gray86')
codeb[,2]=c('white','white','white','white','white',
            'white','white','white','black','black',
            'black','white','white','white','white',
            'black','black','black','white','white',
            'black','white','black','black','black',
            'black','black','black')
b4eq=1:28
nb=length(b4eq)  # nb biomes

load('C5_B4_RCP26.Rdata')
XP=as.numeric(unique(dimnames(B4out)[[2]]))
numax=18 # max number of models
nslice=length(XP)
m=27   # sorties biome4
mt=41  # sorties SB4
ngs=500  # max number of gridcells randomly taken / simulation / scenario
SB4=matrix(NA,numax*nslice*ngs*2,mt)
colnames(SB4)=c("RCP","Sim","Slice","Gridpt","Biome","Count","CO2",
                rownames(F26),rownames(Ins)[iins],"POP","Volc","Solar","domPFT","AET",
                'MTCO','MTWA','E/PE','P-E','GDD5','TANN','PANN','PCOQ','PWAQ','NPPtot',PFTNAMES)
rownames(SB4)=1:(numax*nslice*ngs*2)
layout(matrix(1:4,2,2))
op=par(mai=c(0.8, 0.7, 0.6, 0.2) )

# computation present climate (Mtp): "MTCO0","MTWA0","PCOQ0","PWAQ0"
# we remove the gridpoints with MTWA<0 in Mtp and coord
imtwa=which(apply(Mtp[,1:12],1,max)>0)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
ng=nrow(Mtp)
Z=B4out[1,1,imtwa,]
C0=matrix(NA,ng,mt)
C0[,c(5,17:mt)]=Z[,c(2,4,5:14,3,15:27)]
colnames(C0)=colnames(SB4)
rownames(C0)=1:ng
C0[,"E/PE"]=pmin(100,pmax(0,C0[,"E/PE"]))   # E/PE is delimited by 0 and 100


# RCP2.6
load('C5_B4_RCP26.Rdata')
cat('RCP2.6','\n')
Z=apply(B4out[,1,imtwa,8],2,mean,na.rm=T)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP2.6')
kk=0
num=length(GCM)
print(c(num,nslice))
for (ii in 1:num) {
  for (j in 1:nslice) {
    Z=B4out[ii,j,imtwa,]
    Z[Z[,8]>100,8]=NA
    Z[Z[,8]<0,8]=NA
    zz <- b4eq[Z[,2]]
    h=tabulate(zz)
    if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
    ll=round(pmin(h/ng*ngs,ngs/10),0)
    ll[ll==0] <- NA
    ll=ll+2
    ll[is.na(ll)] <- 0
    for (l in b4eq) {
      if(ll[l]>0) {
        sel=sample(which(zz==l),ll[l])
        for (i in 1:ll[l]) {
          if (Z[sel[i],5]>0) {
          kk=kk+1
          SB4[kk,1]=2.6
          SB4[kk,2]=ii
          SB4[kk,3]=XP[j]
          SB4[kk,4]=sel[i]
          SB4[kk,5]=b4eq[l]
          SB4[kk,6]=h[l]
          SB4[kk,7]=Z[1,1]
          SB4[kk,8:10]=F26[,j+1]
          SB4[kk,11:13]=Ins[iins,7]
          SB4[kk,14]=pop[1,j+1]
          SB4[kk,15:16]=Fpal[2,5:6]   # =midhol =present volc, solar
          SB4[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
          }
        }
      }
    }
  }
}
SB4=SB4[1:kk,]
ir=SB4[SB4[,3]==2006,4]
Z=SB4[SB4[,3]==2006,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP2.6')
print(mean(SB4[,24]))

# RCP4.5
cat('RCP4.5','\n')
load('C5_B4_RCP45.Rdata')
num=length(GCM)
print(c(num,nslice))
SB4b=matrix(NA,num*nslice*ngs*2,mt)
rownames(SB4b)=1:(num*nslice*ngs*2)+nrow(SB4)
Z=apply(B4out[,1,imtwa,8],2,mean,na.rm=T)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP4.5')
kk=0
num=length(GCM)
for (ii in 1:num) {
  for (j in 1:nslice) {
    Z=B4out[ii,j,imtwa,]
    Z[Z[,8]>100,8]=NA
    Z[Z[,8]<0,8]=NA
    zz <- b4eq[Z[,2]]
    h=tabulate(zz)
    if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
    ll=round(pmin(h/ng*ngs,ngs/10),0)
    ll[ll==0] <- NA
    ll=ll+2
    ll[is.na(ll)] <- 0
    for (l in b4eq) {
      if(ll[l]>0) {
        sel=sample(which(zz==l),ll[l])
        for (i in 1:ll[l]) {
          if (Z[sel[i],5]>0) {
          kk=kk+1
          SB4b[kk,1]=4.5
          SB4b[kk,2]=ii
          SB4b[kk,3]=XP[j]
          SB4b[kk,4]=sel[i]
          SB4b[kk,5]=b4eq[l]
          SB4b[kk,6]=h[l]
          SB4b[kk,7]=B4out[ii,j,1,1]
          SB4b[kk,8:10]=F45[,j+1]
          SB4b[kk,11:13]=Ins[iins,7]
          SB4b[kk,14]=pop[2,j+1]
          SB4b[kk,15:16]=Fpal[2,5:6]
          SB4b[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
          }
        }
      }
    }
  }
}
ir=SB4b[SB4b[1:kk,3]==2006,4]
Z=SB4b[SB4b[,3]==2006,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP4.5')
SB4=rbind(SB4,SB4b[1:kk,])
print(mean(SB4b[,24],na.rm=T))

# RCP8.5
cat('RCP8.5','\n')
load('C5_B4_RCP85.Rdata')
num=length(GCM)
print(c(num,nslice))
SB4b=matrix(NA,num*nslice*ngs*2,mt)
rownames(SB4b)=1:(num*nslice*ngs*2)+nrow(SB4)
Z=apply(B4out[,1,imtwa,8],2,mean,na.rm=T)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP8.5')
kk=0
for (ii in 1:num) {
  for (j in 1:nslice) {
    Z=B4out[ii,j,imtwa,]
    Z[Z[,8]>100,8]=NA
    Z[Z[,8]<0,8]=NA
    zz <- b4eq[Z[,2]]
    h=tabulate(zz)
    if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
    ll=round(pmin(h/ng*ngs,ngs/10),0)
    ll[ll==0] <- NA
    ll=ll+2
    ll[is.na(ll)] <- 0
    for (l in b4eq) {
      if(ll[l]>0) {
        sel=sample(which(zz==l),ll[l])
        for (i in 1:ll[l]) {
          if (Z[sel[i],5]>0) {
          kk=kk+1
          SB4b[kk,1]=8.5
          SB4b[kk,2]=ii
          SB4b[kk,3]=XP[j]
          SB4b[kk,4]=sel[i]
          SB4b[kk,5]=b4eq[l]
          SB4b[kk,6]=h[l]
          SB4b[kk,7]=B4out[ii,j,1,1]
          SB4b[kk,8:10]=F85[,j+1]
          SB4b[kk,11:13]=Ins[iins,7]
          SB4b[kk,14]=pop[3,j+1]
          SB4b[kk,15:16]=Fpal[2,5:6]
          SB4b[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
          }
        }
      }
    }
  }
}
SB4=rbind(SB4,SB4b[1:kk,])
ir=SB4b[SB4b[1:kk,3]==2006,4]
Z=SB4b[SB4b[,3]==2006,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / RCP8.5')
print(mean(SB4b[,24],na.rm=T))
rm(SB4b)


# PMIP3 midHolocene
load('P3_B4_midHol.Rdata')
cat('midHol','\n')
num=length(GCM)
print(c(num))
SB4b=matrix(NA,num*ngs*3,mt)
rownames(SB4b)=1:(num*ngs*3)+nrow(SB4)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],B4out[1,imtwa,8],limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / MidHol')
kk=0
for (ii in 1:num) {
  Z=B4out[ii,imtwa,]
  Z[Z[,8]>100,8]=NA   #  E/PE
  Z[Z[,8]<0,8]=NA
  zz <- b4eq[Z[,2]]
  h=tabulate(zz)
  if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
  ll=round(pmin(h/ng*ngs*2,ngs/10),0)
  ll[ll==0] <- NA
  ll=ll+2
  ll[is.na(ll)] <- 0
  for (l in b4eq) {
    if(ll[l]>0) {
      sel=sample(which(zz==l),ll[l])
      for (i in 1:ll[l]) {
        if (Z[sel[i],5]>0) {
        kk=kk+1
        SB4b[kk,1]=-6
        SB4b[kk,2]=ii
        SB4b[kk,3]=-6000
        SB4b[kk,4]=sel[i]
        SB4b[kk,5]=b4eq[l]
        SB4b[kk,6]=h[l]
        SB4b[kk,7]=B4out[ii,1,1]
        SB4b[kk,8:10]=Fpal[2,1:3]
        SB4b[kk,11:13]=Ins[iins,2]
        SB4b[kk,14:16]=Fpal[2,4:6]
        SB4b[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
        }
      }
    }
  }
}
ir=SB4b[1:kk,4]
Z=SB4b[1:kk,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / MidHol')
SB4=rbind(SB4,SB4b[1:kk,])
print(mean(SB4b[1:kk,24],na.rm=T))
rm(SB4b)

# PMIP3 lgm
load('P3_B4_lgm.Rdata')
cat('lgm','\n')
num=length(GCM)
print(c(num))
SB4b=matrix(NA,num*ngs*3,mt)
rownames(SB4b)=1:(num*ngs*3)+nrow(SB4)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],B4out[1,imtwa,8],limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / LGM')
kk=0
for (ii in 1:num) {
  Z=B4out[ii,imtwa,]
  Z[Z[,8]>100,8]=NA
  Z[Z[,8]<0,8]=NA
  zz <- b4eq[Z[,2]]
  h=tabulate(zz)
  if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
  ll=round(pmin(h/ng*ngs*2,ngs/10),0)
  ll[ll==0] <- NA
  ll=ll+2
  ll[is.na(ll)] <- 0
  for (l in b4eq) {
    if(ll[l]>0) {
      sel=sample(which(zz==l),ll[l])
      for (i in 1:ll[l]) {
        if (Z[sel[i],5]>0) {
        kk=kk+1
        SB4b[kk,1]=-21
        SB4b[kk,2]=ii
        SB4b[kk,3]=-21000
        SB4b[kk,4]=sel[i]
        SB4b[kk,5]=b4eq[l]
        SB4b[kk,6]=h[l]
        SB4b[kk,7]=B4out[ii,1,1]
        SB4b[kk,8:10]=Fpal[1,1:3]
        SB4b[kk,11:13]=Ins[iins,1]
        SB4b[kk,14:16]=Fpal[1,4:6]
        SB4b[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
        }
      }
    }
  }
}
SB4=rbind(SB4,SB4b[1:kk,])
ir=SB4b[1:kk,4]
Z=SB4b[1:kk,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / LGM')
print(mean(SB4b[1:kk,24]))
rm(SB4b)

# PMIP3 1K
load('P3_B4_1K.Rdata')
cat('1K','\n')
num=length(GCM)
slice1K=c(1200,1255,1380,1720,1810)
nsl1k=length(slice1K)
print(c(num,nsl1k))
SB4b=matrix(NA,num*nsl1k*ngs*3,mt)
rownames(SB4b)=1:(num*nsl1k*ngs*3)+nrow(SB4)
Mtp=Mtp[imtwa,]
coord=coord[imtwa,]
cartosurf(coord[,1],coord[,2],B4out[1,1,imtwa,8],limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / 1K 1200')
kk=0
for (ii in 1:num) {
  for (is in 1:nsl1k) {
    Z=B4out[ii,is,imtwa,]
    Z[Z[,8]>100,8]=NA
    Z[Z[,8]<0,8]=NA
    zz <- b4eq[Z[,2]]
    h=tabulate(zz)
    if(length(h)<nb) {h=c(h,rep(0,nb-length(h)))}
    ll=round(pmin(h/ng*ngs*2,ngs/10),0)
    ll[ll==0] <- NA
    ll=ll+2
    ll[is.na(ll)] <- 0
    for (l in b4eq) {
      if(ll[l]>0) {
        sel=sample(which(zz==l),ll[l])
        for (i in 1:ll[l]) {
          if (Z[sel[i],5]>0) {
            kk=kk+1
            SB4b[kk,1]=is
            SB4b[kk,2]=ii
            SB4b[kk,3]=slice1K[is]
            SB4b[kk,4]=sel[i]
            SB4b[kk,5]=b4eq[l]
            SB4b[kk,6]=h[l]
            SB4b[kk,7]=B4out[ii,is,1,1]
            SB4b[kk,8:10]=Fpal[2+is,1:3]
            SB4b[kk,11:13]=Ins[iins,2+is]
            SB4b[kk,14:16]=Fpal[2+is,4:6]
            SB4b[kk,17:mt]=Z[sel[i],c(4:14,3,15:m)]
          }
        }
      }
    }
  }
}
SB4=rbind(SB4,SB4b[1:kk,])
ir=SB4b[SB4b[1:kk,3]==1257,4]
Z=SB4b[SB4b[1:kk,3]==1257,21]
cartosurf(coord[ir,1],coord[ir,2],Z,limits=c(-10,50,30,70),zlim=c(0,100)
          ,option="p",coldir=-1,cex=0.25,name='E/PE / 1200')
rm(SB4b)

# correction des outliers: limit to [0,5000]
SB4=na.omit(SB4)
SB4[,"E/PE"]=pmin(100,pmax(0,SB4[,"E/PE"]))  # delimited by [0,100]
SB4[,"P-E"]= pmin(5000,pmax(0,SB4[,"P-E"]))  # delimited by [0,5000]
SB4[,"PANN"]=pmin(5000,pmax(0,SB4[,"PANN"])) # delimited by [0,5000]
# precipitation of the colder quarter
SB4[,"PCOQ"]=pmin(5000,pmax(0,SB4[,"PCOQ"])) # delimited by [0,5000]
# precipitation of the warmer quarter
SB4[,"PWAQ"]=pmin(5000,pmax(0,SB4[,"PWAQ"])) # delimited by [0,5000]
SB4[,"NPPtot"]=pmin(6000,pmax(0,SB4[,"NPPtot"])) # delimited by [0,6000]

# adding random to global variables
dSB4=apply(SB4,2,sd,na.rm=T)
dSB4[14]=dSB4[14]/1000
dSB4[15]=dSB4[15]*10
for (j in 8:16) {
  SB4[,j]=SB4[,j]+0.01*rnorm(nrow(SB4),0,dSB4[j])
}




# aggregation of the pfts
# SB4
Z=matrix(NA,nrow(SB4),8)
colnames(Z)=APFTNAMES
Z[,1]=apply(SB4[,c(29:30)],1,mean,na.rm=T)
Z[,2]=SB4[,31]
Z[,3]=apply(SB4[,c(32:33)],1,mean,na.rm=T)
Z[,4]=apply(SB4[,c(34:35)],1,mean,na.rm=T)
Z[,5]=SB4[,36]
Z[,6]=SB4[,37]
Z[,7]=SB4[,38]
Z[,8]=apply(SB4[,c(39:41)],1,mean,na.rm=T)
SB4=cbind(SB4[,c(1:25,28)],Z)
mt=ncol(SB4)
# C0
Z=matrix(NA,nrow(C0),8)
colnames(Z)=APFTNAMES
Z[,1]=apply(C0[,c(29:30)],1,mean,na.rm=T)
Z[,2]=C0[,31]
Z[,3]=apply(C0[,c(32:33)],1,mean,na.rm=T)
Z[,4]=apply(C0[,c(34:35)],1,mean,na.rm=T)
Z[,5]=C0[,36]
Z[,6]=C0[,37]
Z[,7]=C0[,38]
Z[,8]=apply(C0[,c(39:41)],1,mean,na.rm=T)
C0=cbind(C0[,c(1:25,28)],Z)

# correlation between forcings and bioclimates 
capture.output(print(round(cor(SB4[,18:34],SB4[,8:16]),2),file=RDtxt,append=T))



# calibration GWR: 
# ****************************************************
capture.output(date(),file=RDtxt,append=T)
irg=8:16
idv=c(26:mt,18:25)
il=SB4[,4]
cc=coord[il,]
isel=which(cc[,1]> win[1] & cc[,1]<win[2] & cc[,2]> win[3] & cc[,2]<win[4])
isel1=which(coord[,1]> win[1] & coord[,1]<win[2] & coord[,2]> win[3] & coord[,2]<win[4])
ccw=cc[isel,1:2]
coordw=coord[isel1,]

# PCA of global forcings
X=cbind(SB4[isel,irg])
pcX=prcomp(X,scale.=T)
z=pcX$rotation %*% diag(pcX$sdev)
print('PC loadings')
capture.output(print('PC loadings'),file=RDtxt,append=T)
print(round(z,2))
capture.output(print(round(z,2)),file=RDtxt,append=T)
print('% variance represented by all PCs')
capture.output(print('% variance represented by all PCs'),file=RDtxt,append=T)
print(apply(z^2,1,sum)*100)
capture.output(print(apply(z^2,1,sum)*100),file=RDtxt,append=T)
print('% variance represented by 5 PCs')
capture.output(print('% variance represented by 5 PCs'),file=RDtxt,append=T)
print(round(apply(z[,1:5]^2,1,sum)*100))
capture.output(print(round(apply(z[,1:5]^2,1,sum)*100)),file=RDtxt,append=T)
npcX=5
print(round(apply(z[,1:npcX]^2,1,sum)*100))
capture.output(print(round(apply(z[,1:npcX]^2,1,sum)*100)),file=RDtxt,append=T)
print(paste('% var tot represented by',npcX,'PCs =',round(100*sum(pcX$sdev[1:npcX]^2)/sum(pcX$sdev^2),1)))
capture.output(print(paste('% var tot represented by',npcX,'PCs =',round(100*sum(pcX$sdev[1:npcX]^2)/sum(pcX$sdev^2),1))),file=RDtxt,append=T)

# PFT NPP as function of forcings / GWR
#layout(matrix(1:1,1,1))
#op=par(mai=c(0.6,0.6,0.5, 0.2) )
iosr=c(1:(npcX+1),npcX+2:3,(2*npcX+8):(3*npcX+9))
osr=c("Intercept","PC1","PC2","PC3","PC4","PC5","y","yhat","Intercept_TV","PC1_TV","PC2_TV","PC3_TV","PC4_TV","PC5_TV","Local_R2")
LocCF=array(NA,c(length(isel1),length(iosr),length(idv)))
dimnames(LocCF)[[1]]=rownames(coord)[isel1]
dimnames(LocCF)[[2]]=osr 
dimnames(LocCF)[[3]]=colnames(SB4)[idv]

# on plotte NPP (26), TANN (27), PANN (28)
# for (k in 1:length(idv)) {
for (k in c(1,16,17)) {
  vdn=colnames(SB4)[idv[k]]
  cat('\n')
  cat('Dependent variable',vdn,'\n')
  cat('****************************************','\n')
  cat('\n',file=RDtxt,append=T)
  cat('Dependent variable',vdn,'\n',file=RDtxt,append=T)
  cat('****************************************','\n',file=RDtxt,append=T)
  y0=C0[il,idv[k]]
  Y=SB4[isel,idv[k]]
  data=data.frame(y=Y-y0[isel],pcX$x[,1:npcX],ccw)
  coordinates(data) <- c("long", "lati")
  sr1= gwr.basic (y ~ PC1+PC2+PC3+PC4+PC5, data=data, bw=8)
  print(sr1)
  capture.output(print(sr1),file=RDtxt,append=T)

# plot spatial distribution of the coefficients
  par( mfrow=c(3,3),mai=c(0.3,0.3,0.5, 0.2))
  for (j in 1:npcX) {
    quilt.plot(ccw[,1],ccw[,2],sr1$SDF[[2*npcX+8+j]],main=paste(vdn,'/',colnames(pcX$x)[j],':t-val')
               ,cex.axis=1.2,cex.main=1.6)
    world(add=T)
  }
  quilt.plot(ccw[,1],ccw[,2],sr1$SDF[[(3*npcX+9)]],zlim=c(0,1),main=paste(vdn,'/ Local R2')
             ,cex.axis=1.2,cex.main=1.6)
  world(add=T)
  
  zlim=range(sr1$SDF$yhat)
  quilt.plot(ccw[X[,1]<290,1],ccw[X[,1]<290,2],sr1$SDF$yhat[X[,1]<290],zlim=zlim,main=paste(vdn,'- Est CO2<290')
             ,cex.axis=1.2,cex.main=1.6)
  world(add=T)
  print(mean(sr1$SDF$yhat[X[,1]<290]))
  capture.output(cat('- Est CO2<290',mean(sr1$SDF$yhat[X[,1]<290])),file=RDtxt,append=T)
  quilt.plot(ccw[X[,1]>600,1],ccw[X[,1]>600,2],sr1$SDF$yhat[X[,1]>600],zlim=zlim,main=paste(vdn,'- Est CO2>600')
             ,cex.axis=1.2,cex.main=1.6)
  world(add=T)
  print(mean(sr1$SDF$yhat[X[,1]>600]))
  capture.output(cat(mean(sr1$SDF$yhat[X[,1]>600])),file=RDtxt,append=T)
  
  #quilt.plot(ccw[X[,8]<0.2,1],ccw[X[,8]<0.2,2],sr1$SDF$yhat[X[,8]<0.2],main=paste(vdn,'- Est Volc<0.2'))
  #world(add=T)
  print(mean(sr1$SDF$yhat[X[,8]<0.2]))
  capture.output(cat(mean(sr1$SDF$yhat[X[,8]<0.2])),file=RDtxt,append=T)
  #quilt.plot(ccw[X[,8]>0.4,1],ccw[X[,8]>0.4,2],sr1$SDF$yhat[X[,8]>0.4],main=paste(vdn,'- Est Volc>0.4'))
  #world(add=T)
  print(mean(sr1$SDF$yhat[X[,8]>0.4]))
  capture.output(cat(mean(sr1$SDF$yhat[X[,8]>0.4])),file=RDtxt,append=T)
  
  #quilt.plot(ccw[X[,9]<200,1],ccw[X[,9]<200,2],sr1$SDF$yhat[X[,9]<200],main=paste(vdn,'- Est Solar<200'))
  #world(add=T)
  print(mean(sr1$SDF$yhat[X[,9]<200]))
  capture.output(cat(mean(sr1$SDF$yhat[X[,9]<200])),file=RDtxt,append=T)
  #quilt.plot(ccw[X[,9]>400,1],ccw[X[,9]>400,2],sr1$SDF$yhat[X[,9]>400],main=paste(vdn,'- Est Solar>400'))
  #world(add=T)
  print(mean(sr1$SDF$yhat[X[,9]>400]))
  capture.output(cat(mean(sr1$SDF$yhat[X[,9]>400])),file=RDtxt,append=T)
  
# save spatial distribution of the coefficients
  for (j in 1:length(iosr)) {
    zg=interpp(ccw[,1],ccw[,2],sr1$SDF[[iosr[j]]],duplicate='mean',coordw[,1],coordw[,2])
    LocCF[,j,k]=zg$z
  }
  
#  residuals analysis for four CO2 situations
  xlim=range(sr1$SDF$y)
  plot(xlim,xlim,cex=0.1,xlab='Obs',ylab='Pred',main=paste(vdn,'/ Pred vs Obs'),cex.axis=1.2,cex.main=1.6)
  points(sr1$SDF$y[X[,1]<250],sr1$SDF$yhat[X[,1]<250],cex=0.2,col='blue')
  points(sr1$SDF$y[X[,1]>=250 & X[,1]<400],sr1$SDF$yhat[X[,1]>=250 & X[,1]<400],cex=0.2,col='green')
  points(sr1$SDF$y[X[,1]>=400 & X[,1]<600],sr1$SDF$yhat[X[,1]>=400 & X[,1]<600],cex=0.2,col='orange')
  points(sr1$SDF$y[X[,1]>=600],sr1$SDF$yhat[X[,1]>=600],cex=0.2,col='red')
  abline(a=0,b=1,col='black',lwd=1)  
  text(xlim[1]+(xlim[2]-xlim[1])/10,xlim[2]-(xlim[2]-xlim[1])/10,'colors=CO2')

    
# final save  
#  save.image(RDfile)
}

#save.image(RDfile)
save(Mtp,coord,LocCF,C0,coordw,SB4,file=RDfil2)
capture.output(date(),file=RDtxt,append=T)
dev.off()
