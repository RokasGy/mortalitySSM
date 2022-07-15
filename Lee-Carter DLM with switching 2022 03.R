#################### Libraries 

library(dlm)

#################### Data import

#Set age and year vectors
Age<-seq(from=27.5, by =5, length=10) #10 age groups covering 25-74 years
AgeGroup<-c("25-29","30-34","35-39","40-44","45-49","50-54", 
"55-59","60-64","65-69","70-74")
Year.SEA<-1900:2017; per.SEA<-118  #118 years

#Set working directory for input/output
#SPECIFY YOUR WORKING DIRECTORY HERE
setwd("C:\\....\\")

#THE MODEL USES (SWEDISH) MORTALITY DATA DOWNLOADED FROM MORTALITY.ORG
#Data file used in the calculations contained data up to 2018
E.SEA<- read.table(".\\Data\\Exposures_5x1.txt",
skip=2, header=TRUE,colClasses = "character")
E.SEA<-E.SEA[,c(1,2,5)]
E.SEA$Year<-as.numeric(E.SEA$Year)
E.SEA$Age<-as.factor(E.SEA$Age)
E.SEA$Total<-as.numeric(E.SEA$Total)
E.SEA<-E.SEA[E.SEA$Age %in% AgeGroup,]
E.SEA<-E.SEA[E.SEA$Year>1899&E.SEA$Year<2018,]
EM.SEA<-matrix(E.SEA$Total,nrow=10,dimnames=list(AgeGroup,Year.SEA))

D.SEA<- read.table(".\\Data\\Deaths_5x1.txt",
skip=2, header=TRUE,colClasses = "character")
D.SEA<-D.SEA[,c(1,2,5)]
D.SEA$Year<-as.numeric(D.SEA$Year)
D.SEA$Age<-as.factor(D.SEA$Age)
D.SEA$Total<-round(as.numeric(D.SEA$Total),0)
D.SEA<-D.SEA[D.SEA$Age %in% AgeGroup,]
D.SEA<-D.SEA[D.SEA$Year>1899&D.SEA$Year<2018,]
DM.SEA<-matrix(D.SEA$Total,nrow=10,dimnames=list(AgeGroup,Year.SEA))

Mx.SEA<-D.SEA[,c(1,2)]
Mx.SEA$Mx<-D.SEA$Total/E.SEA$Total
MxM.SEA<-matrix(Mx.SEA$Mx,nrow=10,dimnames=list(AgeGroup,Year.SEA)) #matrix of forces of mortality






#################### Classical LC

a.LC.SEA <- rowMeans(log(MxM.SEA))
MxMc.SEA<-MxM.SEA
	for(i in 1:10) MxMc.SEA[i,] <- log(MxM.SEA[i,]) - a.LC.SEA[i]  #centered MxM
LC.SEA<- svd(MxMc.SEA, 1, 1)
b.LC.SEA <- LC.SEA$u/sum(LC.SEA$u) #normalisation of parameters
k.init.SEA <- LC.SEA$v*sum(LC.SEA$u)*LC.SEA$d[1]

k.SEA<-ts(k.init.SEA, start=1900,freq=1)
ts.SEA<-arima(k.SEA, order=c(0,1,0),xreg=1:per.SEA)
drift.SEA<-as.numeric(ts.SEA$coef)
s2.SEA<-ts.SEA$sigma2 
s2.drift.SEA<-as.numeric(ts.SEA$var.coef)

k.SEA.fit<-predict(ts.SEA,newxreg=1:per.SEA)
Mx.SEA.fit<-exp(b.LC.SEA%*%t(k.SEA.fit$pred)+a.LC.SEA)
res<-log(MxM.SEA)-log(Mx.SEA.fit)
s2.V.SEA<-sum(res^2)/(per.SEA*10)


#### Forecasting and simulation

k.SEA.pred<-predict(ts.SEA,n.ahead=30,newxreg=(per.SEA+1):(per.SEA+30)) 
Mx.SEA.f<-exp(b.LC.SEA%*%t(k.SEA.pred$pred)+a.LC.SEA)
colnames(Mx.SEA.f)<-2018:2047
#Qx.SEA.f<-1-exp(-Mx.SEA.f)

# k(t) simulation for 30 years
set.seed(101)
k.SEA.sim = matrix(nrow=20000, ncol=30)
Mx.sim<-matrix(nrow=10,ncol=30)
Mx.SEA.sim<-matrix(nrow=20000, ncol=300)

for (i in 1:20000){
	m<-rnorm(n=1,mean=drift.SEA, sd=sqrt(s2.drift.SEA))
	x<-k.SEA[per.SEA]+cumsum(rnorm(n=30, mean=m, sd = sqrt(s2.SEA))) 
	k.SEA.sim[i,]<-x

mu=a.LC.SEA+b.LC.SEA%*%t(x)
	for (m in 1:30){
	Mx.sim[,m]<-mu[,m]+rnorm(n=10,mean=0,sd=sqrt(s2.V.SEA))
	}
Mx.SEA.sim[i,]<-as.vector(exp(Mx.sim))
}

CI95.k.SEA.sim<-apply(k.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.k.SEA.sim<-apply(k.SEA.sim,2,quantile, probs = c(.005,.5,.995))

CI95.Mx.SEA.sim<-apply(Mx.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.Mx.SEA.sim<-apply(Mx.SEA.sim,2,quantile, probs = c(.005,.5,.995))






#################### DLM LC

buildFunB1 <- function(x, start.x, start.V) {
A <- matrix(c(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]), 
	ncol = 1) #matrix of factors loadings
FF <-  A %*% matrix(c(0,0,1),nrow=1)
GG <- matrix(c(1,0,1,0,1,100,0,0,1),ncol=3)
JGG<-matrix(c(0,0,0,0,0,1,0,0,0), ncol=3)
W <- diag(c(0,0,exp(x[11])), ncol = 3)
JW<-diag(c(0,0,2), ncol = 3)
X<-matrix(rep(0,per.SEA),ncol=1)
V<-diag(rep(exp(x[12]),10),ncol=10) 
dlmSEA<-dlm(FF=FF,V=V,GG=GG,JGG=JGG,W=W,X=X,m0=c(x[13],100,start.x),C0=diag(c(0,0,start.V))) 
return(dlmSEA)
}

chain<-5 #number of chains
mcmc<-5000 #number of itterations
burn<-1000
every<-1

x<-rep(0.1,14)
x0<-matrix(nrow=chain,ncol=14)
x0[1,]<-rep(0.1,14)	#Par 14 not used in SEA case 
x0[2,]<-c(rep(0.1,10),1,0.1,0.1,0.1)
x0[3,]<-c(rep(0.1,10),0.1,1,0.1,0.1)
x0[4,]<-c(rep(0.1,10),0.1,0.1,1,0.1)
x0[5,]<-c(rep(0.1,10),0.5,0.5,0.5,0.1)
start.x<-9  #starting parameter k value 
start.V<-10

xx1.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=14),simplify=FALSE)
kk1.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=(per.SEA+1)),simplify=FALSE)

prior.m<-0
prior.m.prec<-1/5
prior.g1<-2.1
prior.g2<-0.1   
prior.b<-0.1 #equal for each age group 
prior.b.prec<-1/5
 
set.seed(101)
for (ii in 1:chain) {
	x<-x0[ii,]
	i.save<-0

   for (i in 1:mcmc) { 
	
	# 1. Generate sample of states
	mod=buildFunB1(x, start.x, start.V)
	modFilt <- dlmFilter(t(as.matrix(MxMc.SEA)), mod, simplify=TRUE)
	k <- dlmBSample(modFilt)  


	# 2. Sample parameter values

	# 2.1. sample drifts

	xs<-rep(1,per.SEA)/sqrt(exp(x[11]))
	dk<-k[-1,3]-k[-nrow(k),3]
	dks<-dk/sqrt(exp(x[11]))

	mu0<-rnorm(n=1,mean=(prior.m.prec*prior.m+sum(xs*dks))/(prior.m.prec+sum(xs*xs)), 
		sd=sqrt(1/(prior.m.prec+sum(xs*xs))))


	# 2.2. sample W
	res1<-dk-mu0
	w0.0<-1/rgamma(n=1, prior.g1+per.SEA/2, prior.g2+0.5*sum(res1^2))
	

	# 2.3. sample betas
	kk<-t(k[-1,3])%*%k[-1,3]
	b0<-1:10
	for(j in 1:10){	
		b0[j]<-rnorm(n=1, mean=(prior.b.prec*prior.b*exp(x[12])+(t(k[-1,3])%*%MxMc.SEA[j,]))/
			(prior.b.prec*exp(x[12])+kk),	sd=sqrt(exp(x[12])/(exp(x[12])*prior.b.prec+kk)))      
	}


	# 2.4. sample V
	res3 <- MxMc.SEA-b0%*%t(k[-1,3])
	v<-1/rgamma(n=1, prior.g1+per.SEA*10/2, prior.g2+0.5*sum(res3^2))    


	# 3.6. reweight betas and k(t) parameters

	bb<-sum(b0)
	b<-b0/bb
	mu<-mu0*bb
	w.0<-w0.0*(bb^2)
	
		
	# 4. Populate updated paramters to x
	x<-c(b,log(w.0),log(v),mu,100)    

	# 5. Save selected itterations

	if(i>burn && !(i %% every)){
		i.save<-i.save+1
		kk1.SEA[[ii]][i.save,]<-k[,3]
		xx1.SEA[[ii]][i.save,]<-x
	}
   }
}


#### Parameter summary

par.b1.m<-matrix(nrow=chain*2, ncol=14)
par.b1.var<-matrix(nrow=chain*2, ncol=14)

for (jj in 1:chain) {
	xx<-xx1.SEA[[jj]]
	xx[,11:12]<-exp(xx[,11:12])
	par.b1.m[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,mean)
	par.b1.m[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,mean)
	par.b1.var[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,var)
	par.b1.var[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,var)
}
#examine for significant deviations between results of different half-chains

par.b1<-apply(par.b1.m,2,mean)


#### MCMC diagnostics

#Calculation of B and W statistics
B.mc1<-rep(0,14)
W.mc1<-rep(0,14)

n=(mcmc-burn)/2
m=chain*2
for (ii in 1:14) {
	B.mc1[ii]<-sum((par.b1.m[,ii]-par.b1[ii])^2)*n/(m-1)  #var between chains
	W.mc1[ii]<-sum(par.b1.var[,ii])/m  #within chain variance
}

V.mc1<-(n-1)*W.mc1/n+B.mc1/n
sqrt(V.mc1)
(Rhat<-sqrt(V.mc1/W.mc1))  #R hat statistics calculation

xx1.SEA.m<-do.call(rbind,xx1.SEA)
acf<-apply(xx1.SEA.m,2,function(x) acf(x,plot=FALSE)$acf[-1])

#vector of indexes of acfs 
acf.d<-acf[-nrow(acf),]+acf[-1,]
ind<-apply(acf.d,2, function(x) ifelse(is.na(which(x<=0)[1]),nrow(acf),which(x<=0)[1]))

Eff.mc1<-rep(0,14)

for (ii in 1:14) {
	x2<-acf[1:ind[ii],ii]
	Eff.mc1[ii]<-m*n/(1+2*sum(x2))  
	}

Eff.mc1


#### Bayesian forecasting based on posterior samples

nsim<-20000
Mx.SEA2.sim<-matrix(nrow=nsim,ncol=300)
colnames(Mx.SEA2.sim)<-paste(rep(1:30, each=10),AgeGroup)
k.SEA2.sim<-matrix(nrow=nsim,ncol=30)
nr<-nrow(xx1.SEA.m)
kk1.SEA.m<-do.call(rbind,kk1.SEA)
start.x<-colSums(kk1.SEA.m)[length(colSums(kk1.SEA.m))]/nr
start.V<-0.0001
k.sim<-c(start.x,rep(100,30))
Mx.sim<-matrix(nrow=10,ncol=30)

set.seed(105)
for(j in 1:nsim){
samp.ind<-sample(x=nr,1,rep=TRUE)
	for (k in 1:30){
	k.sim[k+1]<-k.sim[k]+rnorm(1,mean=xx1.SEA.m[samp.ind,13],
		sd=sqrt(exp(xx1.SEA.m[samp.ind,11])))
	}
k.SEA2.sim[j,]<-k.sim[-1]
mu=a.LC.SEA+xx1.SEA.m[samp.ind,1:10]%*%t(k.sim[-1])
	for (m in 1:30){
	Mx.sim[,m]<-mu[,m]+rnorm(n=10,mean=0,sd=sqrt(exp(xx1.SEA.m[samp.ind,12])))
	}
Mx.SEA2.sim[j,]<-as.vector(exp(Mx.sim))
}

CI95.k.SEA2.sim<-apply(k.SEA2.sim,2,quantile, probs = c(.025,.5,.975))
CI99.k.SEA2.sim<-apply(k.SEA2.sim,2,quantile, probs = c(.005,.5,.995))

CI95.Mx.SEA2.sim<-apply(Mx.SEA2.sim,2,quantile, probs = c(.025,.5,.975))
CI99.Mx.SEA2.sim<-apply(Mx.SEA2.sim,2,quantile, probs = c(.005,.5,.995))






#################### DLM LC with switching k(t) variance

buildFunB2 <- function(x, vecX, start.x, start.V) {
A <- matrix(c(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]), 
	ncol = 1) #matrix of factors loadings
FF <-  A %*% matrix(c(0,0,1),nrow=1)
GG <- matrix(c(1,0,1,0,1,100,0,0,1),ncol=3)
JGG<-matrix(c(0,0,0,0,0,1,0,0,0), ncol=3)
W <- diag(c(0,0,100), ncol = 3)
JW<-diag(c(0,0,2), ncol = 3)
X<-cbind(rep(0,per.SEA),exp(x[12])*matrix(vecX,ncol=1)+exp(x[11])*(1-matrix(vecX,ncol=1)))
V<-diag(rep(exp(x[13]),10),ncol=10) 
dlmSEA<-dlm(FF=FF,V=V,GG=GG,JGG=JGG,W=W,JW=JW,X=X,m0=c(x[14],100,start.x),C0=diag(c(0,0,start.V))) 
return(dlmSEA)
}

chain<-5 
mcmc<-5000
burn<-1000
every<-1

x<-c(rep(0.1,14),100)	# Par 15 not used in SEA case 
x0<-matrix(nrow=chain,ncol=15)
x0[1,]<-c(rep(0.1,14),100)		#Par 15 not used in SEA case 
x0[2,]<-c(rep(0.1,10),0.05,0.1,0.1,0.1,100)
x0[3,]<-c(rep(0.1,10),0.5,1,0.1,0.1,100)
x0[4,]<-c(rep(0.1,10),0.1,0.1,1,0.1,100)
x0[5,]<-c(rep(0.1,10),0.1,0.1,0.1,1,100)
start.x<-9  
start.V<-10

xx2.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=17),simplify=FALSE)
kk2.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=(per.SEA+1)),simplify=FALSE)
Im.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=per.SEA),simplify=FALSE)
tbl.SEA<-replicate((mcmc-burn)*chain/every, matrix(nrow=2, ncol=2),simplify=FALSE)

p<-matrix(nrow=2, ncol=(per.SEA+1)) 
ind<-matrix(nrow=2, ncol=per.SEA)

prior.beta<-1  
prior.beta1<-1  
prior.beta2<-1 
prior.m<-0
prior.m.prec<-1/5
prior.g1<-2.1
prior.g2<-0.1   
prior.b<-0.1 #equal for each age group 
prior.b.prec<-1/5


set.seed(102)
for (ii in 1:chain) {
	x<-x0[ii,]
	i.save<-0
	i0<-0.5  #transition probability from 0 to 0
	i1<-0.5  #transition probability from 1 to 1
	vecX<-rep(0, per.SEA)

   for (i in 1:mcmc) { 
	
	# 1. Generate sample of states
	mod=buildFunB2(x, vecX,start.x, start.V)
	modFilt <- dlmFilter(t(as.matrix(MxMc.SEA)), mod, simplify=TRUE)
	k <- dlmBSample(modFilt)  


	# 2. Generate mixture index
	# 2.1. generate probabilities
	p[1,1]<-(1-i1)/(2-i1-i0)  
	p[2,1]<-(1-i0)/(2-i1-i0)
	
	for (j in 1:per.SEA) { 
		ind[1,j]<-p[1,j]*i0+p[2,j]*(1-i1)  #prior probabilities of states
		ind[2,j]<-p[2,j]*i1+p[1,j]*(1-i0)  
		
		lk1<-dnorm(x=(k[j+1,3]-k[j,3]),mean=k[j,1],sd=sqrt(exp(x[11])))  
		lk2<-dnorm(x=(k[j+1,3]-k[j,3]),mean=k[j,1],sd=sqrt(exp(x[12])))

		p0<-ind[1,j]*lk1  #posterior probabilities
		p1<-ind[2,j]*lk2

      	p[1,j+1]<-p0/(p0+p1) #normalised posterior probabilities
		p[2,j+1]<-p1/(p0+p1)
		}  

	
	# 2.2. perform conditional backwards sampling
	vecX[per.SEA]<-rbinom(n=1, size=1, prob=p[2,(per.SEA+1)])  #p of getting 1
	for (h in (per.SEA-1):1) {
		v0=i0*(1-vecX[h+1])+(1-i0)*vecX[h+1] #if vecX[h+1]=1, v0 gives p of moving from state 0 to 1
								#if vecX[h+1]=0, v0 gives p of staying in state 0 (1-i0)
		v1=i1*(vecX[h+1])-(1-i1)*(vecX[h+1]-1) #if vecX[h+1]=1, v1 gives p of staying in state 1 (1-i1)
								#if vecX[h+1]=0, v1 gives p of moving from state 1 to 0
		p0=v0*p[1,h+1]   #joint p of being in state 0 before transition
 		p1=v1*p[2,h+1]   #joint p of being in state 1 before transition
		p0.s<-p0/(p0+p1)
		p1.s<-p1/(p0+p1)
		vecX[h]<-rbinom(n=1, size=1, prob=p1.s)
		}
	

	# 3. Sample parameter values

	# 3.1. sample transition probabilities

	tr<-matrix(rep(0,4),2,2)
	tryCatch(
		tr[1,1]<-table(head(vecX,-1),tail(vecX,-1))["0","0"],
		error=function(e) {tr[1,1]<-0})
	tryCatch(
		tr[1,2]<-table(head(vecX,-1),tail(vecX,-1))["0","1"], #transitions from state 1 to 0
		error=function(e) {tr[1,2]<-0})
	tryCatch(
		tr[2,2]<-table(head(vecX,-1),tail(vecX,-1))["1","1"],
		error=function(e) {tr[2,2]<-0})
	tryCatch(
		tr[2,1]<-table(head(vecX,-1),tail(vecX,-1))["1","0"], #transitions from state 0 to 1
		error=function(e) {tr[2,1]<-0})
	
	i0<-rbeta(n=1, prior.beta+tr[1,1], prior.beta1+tr[2,1])
	i1<-rbeta(n=1, prior.beta2+tr[2,2], prior.beta1+tr[1,2])


	# 3.2. sample drifts

	aa<-sqrt(exp(x[12])*matrix(vecX,ncol=1)+exp(x[11])*(1-matrix(vecX,ncol=1))) #vector of sd's 
	xs<-1/aa
	dk<-k[-1,3]-k[-nrow(k),3]
	dks<-dk/aa

	mu0<-rnorm(n=1,mean=(prior.m.prec*prior.m+sum(xs*dks))/(prior.m.prec+sum(xs*xs)), 
		sd=sqrt(1/(prior.m.prec+sum(xs*xs))))


	# 3.3. sample W

	hh=exp(x[12])/exp(x[11])-1    
	dkss<-dk/sqrt(1+hh*vecX)
	res1<-dkss-mu0/sqrt(1+hh*vecX)

	w0.0<-1/rgamma(n=1, prior.g1+per.SEA/2, prior.g2+0.5*sum(res1^2))

	dkss1<-(dk/sqrt(w0.0)-mu0/sqrt(w0.0))*vecX
	res2<-dkss1[dkss1!=0]
	
	hh<-1/rgamma(n=1, prior.g1+length(res2)/2, prior.g2+0.5*sum(res2^2))-1
	
	ifelse(hh>0,
		w1.0<-(1+hh)*w0.0,
		w1.0<-exp(x[12])*w0.0/exp(x[11]))
		

	# 3.4. sample betas

	kk<-t(k[-1,3])%*%k[-1,3]
	b0<-1:10
	for(j in 1:10){	
		b0[j]<-rnorm(n=1, mean=(prior.b.prec*prior.b*exp(x[13])+(t(k[-1,3])%*%MxMc.SEA[j,]))
			/(prior.b.prec*exp(x[13])+kk),sd=sqrt(exp(x[13])/(prior.b.prec*exp(x[13])+kk)))      
	}


	# 3.5. sample V
	res3 <- MxMc.SEA-b0%*%t(k[-1,3])
	v<-1/rgamma(n=1, prior.g1+per.SEA*10/2, prior.g2+0.5*sum(res3^2))    


	# 3.6. reweight betas and k(t) parameters

	bb<-sum(b0)
	b<-b0/bb
	mu<-mu0*bb
	w.0<-w0.0*(bb^2)
	w.1<-w1.0*(bb^2)
	
		
	# 4. Populate updated paramters to x

	x<-c(b,log(w.0),log(w.1),log(v),mu,100)    

	# 5. Save selected itterations

	if(i>burn && !(i %% every)){
		i.save<-i.save+1
		kk2.SEA[[ii]][i.save,]<-k[,3]
		xx2.SEA[[ii]][i.save,]<-c(x,i0,i1)
		Im.SEA[[ii]][i.save,]<-vecX
		tbl.SEA[[(ii-1)*(mcmc-burn)+i.save]]<-tr
	}
   }
}


## Parameter summary

par.b2.m<-matrix(nrow=chain*2, ncol=17)
par.b2.var<-matrix(nrow=chain*2, ncol=17)


for (jj in 1:chain) {
	xx<-xx2.SEA[[jj]]
	xx[,11:13]<-exp(xx[,11:13])
	par.b2.m[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,mean)
	par.b2.m[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,mean)
	par.b2.var[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,var)
	par.b2.var[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,var)
}
#check for significant deviations between results of different half-chains

par.b2<-apply(par.b2.m,2,mean)


## MCMC diagnostics

#Calculation of B and W statistics
B.mc2<-rep(0,17)
W.mc2<-rep(0,17)

n=(mcmc-burn)/2
m=chain*2
for (ii in 1:17) {
	B.mc2[ii]<-sum((par.b2.m[,ii]-par.b2[ii])^2)*n/(m-1)  #var between chains
	W.mc2[ii]<-sum(par.b2.var[,ii])/m  #within chain variance
}

V.mc2<-(n-1)*W.mc2/n+B.mc2/n
sqrt(V.mc2)
(Rhat<-sqrt(V.mc2/W.mc2))  

xx2.SEA.m<-do.call(rbind,xx2.SEA)
acf<-apply(xx2.SEA.m,2,function(x) acf(x,plot=FALSE)$acf[-1])

#vector of indexes of acfs 
acf.d<-acf[-nrow(acf),]+acf[-1,]
ind<-apply(acf.d,2, function(x) ifelse(is.na(which(x<=0)[1]),nrow(acf),which(x<=0)[1]))

Eff.mc2<-rep(0,17)

for (ii in 1:17) {
	x2<-acf[1:ind[ii],ii]
	Eff.mc2[ii]<-m*n/(1+2*sum(x2))  
	}

Eff.mc2



#Bayesian forecasting based on posterior samples

kk2.SEA.m<-do.call(rbind,kk2.SEA)
nsim<-20000
Mx.SEA3.sim<-matrix(nrow=nsim,ncol=300)
colnames(Mx.SEA3.sim)<-paste(rep(1:30, each=10),AgeGroup)
p.sim<-matrix(nrow=2,ncol=31)
vecX.sim<-rep(10,30)
Im.SEA3.sim<-matrix(nrow=nsim,ncol=30)
k.SEA3.sim<-matrix(nrow=nsim,ncol=30)
nr<-nrow(xx2.SEA.m)
start.x<-colSums(kk2.SEA.m)[length(colSums(kk2.SEA.m))]/nr
k.sim<-c(start.x,rep(100,30))
Mx.sim<-matrix(nrow=10,ncol=30)


set.seed(107)
for(j in 1:nsim){
	samp.ind<-sample(x=nr,1,rep=TRUE)
	
	#simulate the regime vector
	i0<-xx2.SEA.m[samp.ind,16]
	i1<-xx2.SEA.m[samp.ind,17]
	p.sim[1,]<-(1-i1)/(2-i1-i0)  
	p.sim[2,]<-(1-i0)/(2-i1-i0)
	vecX.sim[30]<-rbinom(n=1, size=1, prob=p.sim[2,31])  #p of getting 1
	for (h in 29:1) {
		v0=i0*(1-vecX.sim[h+1])+(1-i0)*vecX.sim[h+1] 
		v1=i1*(vecX.sim[h+1])-(1-i1)*(vecX.sim[h+1]-1) 
		p0=v0*p.sim[1,h+1]   
 		p1=v1*p.sim[2,h+1]  
		p0.s<-p0/(p0+p1)
		p1.s<-p1/(p0+p1)
		vecX.sim[h]<-rbinom(n=1, size=1, prob=p1.s)
		}
	Im.SEA3.sim[j,]<-vecX.sim
	
	#simulate k(t)
	var.vect<-vecX.sim*xx2.SEA.m[samp.ind,12]+(1-vecX.sim)*xx2.SEA.m[samp.ind,11]
	for (k in 1:30){
		k.sim[k+1]<-k.sim[k]+rnorm(1,mean=xx2.SEA.m[samp.ind,14],sd=sqrt(exp(var.vect[k])))
		}
	k.SEA3.sim[j,]<-k.sim[-1]

	#simulate mortality rates
	mu=a.LC.SEA+xx2.SEA.m[samp.ind,1:10]%*%t(k.sim[-1])
	for (m in 1:30){
		Mx.sim[,m]<-mu[,m]+rnorm(n=10,mean=0,sd=sqrt(exp(xx2.SEA.m[samp.ind,13])))
		}
	Mx.SEA3.sim[j,]<-as.vector(exp(Mx.sim))
	}

CI95.k.SEA3.sim<-apply(k.SEA3.sim,2,quantile, probs = c(.025,.5,.975))
CI99.k.SEA3.sim<-apply(k.SEA3.sim,2,quantile, probs = c(.005,.5,.995))

CI95.Mx.SEA3.sim<-apply(Mx.SEA3.sim,2,quantile, probs = c(.025,.5,.975))
CI99.Mx.SEA3.sim<-apply(Mx.SEA3.sim,2,quantile, probs = c(.005,.5,.995))


