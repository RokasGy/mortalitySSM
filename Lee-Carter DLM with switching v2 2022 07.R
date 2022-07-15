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
#setwd("C:\\....\\")


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




#################### DLM LC

buildFunB1 <- function(x, start.x, start.V) {
FF <- matrix(c(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],rep(0,20),x[11],x[12],x[13],
	x[14],x[15],x[16],x[17],x[18],x[19],x[20]),ncol = 4)
GG <- matrix(c(1,0,0,0,0,1,0,1,0,0,1,100,0,0,0,1),ncol=4)
JGG<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0), ncol=4)
W <- diag(c(0,0,0,exp(x[21])), ncol = 4)
X<-matrix(rep(0,per.SEA),ncol=1)
V<-diag(rep(exp(x[22]),10),ncol=10) 
dlmSEA<-dlm(FF=FF,V=V,GG=GG,JGG=JGG,W=W,X=X,m0=c(1,x[23],x[24],start.x),
	C0=diag(c(0.000001,0.000001,0.000001,start.V))) 
return(dlmSEA)
}

chain<-5  #number of chains
mcmc<-5000 #number of itterations
burn<-1000
every<-1

x<-c(rep(0.1,23),100)	 #Par 24 not used in SEA case
x0<-matrix(nrow=chain,ncol=24)
x0[1,]<-rep(0.1,24)	 
x0[2,]<-c(rep(0.1,20),1,0.1,0.1,100)
x0[3,]<-c(rep(0.1,20),0.1,1,0.1,100)
x0[4,]<-c(rep(0.1,20),0.1,0.1,1,100)
x0[5,]<-c(rep(0.1,20),0.5,0.5,0.5,100)
start.x<-9  #starting parameter k value
start.V<-10

xx1.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=24),simplify=FALSE)
kk1.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=(per.SEA+1)),simplify=FALSE)

prior.m<-0
prior.m.prec<-1/5
prior.g1<-2.1
prior.g2<-0.1   
prior.b<-c(-4,0.1) #equal for each age group 
prior.b.prec<-matrix(c(1/5,0,0,1/5),ncol=2)
 
set.seed(101)
for (ii in 1:chain) {
	x<-x0[ii,]
	i.save<-0

   for (i in 1:mcmc) { 
	
	# 1. Generate sample of states
	mod=buildFunB1(x, start.x, start.V)
	modFilt <- dlmFilter(t(as.matrix(log(MxM.SEA))), mod, simplify=TRUE)
	k <- dlmBSample(modFilt)  


	# 2. Sample parameter values

	# 2.1. sample drifts

	xs<-rep(1,per.SEA)/sqrt(exp(x[21]))
	dk<-k[-1,4]-k[-nrow(k),4]
	dks<-dk/sqrt(exp(x[21]))

	mu0<-rnorm(n=1,mean=(prior.m.prec*prior.m+sum(xs*dks))/(prior.m.prec+sum(xs*xs)), 
		sd=sqrt(1/(prior.m.prec+sum(xs*xs))))


	# 2.2. sample W
	res1<-dk-mu0
	w0.0<-1/rgamma(n=1, prior.g1+per.SEA/2, prior.g2+0.5*sum(res1^2))


	# 2.3. sample alphas and betas
	kk<-rbind(rep(1,per.SEA),k[-1,4])%*%cbind(rep(1,per.SEA),k[-1,4])
	b0<-cbind(1:10,1:10)
	for(j in 1:10){	
		
		Mu=solve(prior.b.prec+(1/exp(x[22]))*kk)%*%(prior.b.prec%*%prior.b+(1/exp(x[22]))*
			rbind(rep(1,59),k[-1,4])%*%log(MxM.SEA[j,]))
		Sigma=solve(prior.b.prec+(1/exp(x[22]))*kk)

		b01<-rnorm(n=1,mean=Mu[2],sd=sqrt(Sigma[2,2]))
		b00<-rnorm(n=1,mean=Mu[1]+Sigma[1,2]*(b01-Mu[2])/Sigma[2,2],
			sd=sqrt(Sigma[1,1]-Sigma[1,2]*Sigma[2,1]/Sigma[2,2]))
		b0[j,]<-c(b00,b01)
  
	}


	# 2.4. sample V
	res3 <- log(MxM.SEA)-b0[,1]-b0[,2]%*%t(k[-1,4])
	v<-1/rgamma(n=1, prior.g1+per.SEA*10/2, prior.g2+0.5*sum(res3^2))    


	# 3.6. reweight betas and k(t) parameters
	bb<-sum(b0[,2])
	b<-b0[,2]/bb
	a<-b0[,1]+mean(k[-1,4])*b0[,2]
	mu<-(mu0-mean(k[-1,4]))*bb
	w.0<-w0.0*(bb^2)
	
		
	# 4. Populate updated paramters to x
	x<-c(a,b,log(w.0),log(v),mu,100)    

	# 5. Save selected itterations

	if(i>burn && !(i %% every)){
		i.save<-i.save+1
		kk1.SEA[[ii]][i.save,]<-k[,4]
		xx1.SEA[[ii]][i.save,]<-x
	}
   }
}

## Parameter summary

par.b1.m<-matrix(nrow=chain*2, ncol=24)
par.b1.var<-matrix(nrow=chain*2, ncol=24)

for (jj in 1:chain) {
	xx<-xx1.SEA[[jj]]
	xx[,21:22]<-exp(xx[,21:22])
	par.b1.m[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,mean)
	par.b1.m[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,mean)
	par.b1.var[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,var)
	par.b1.var[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,var)
}
#examine for significant deviations between results of different half-chains

par.b1<-apply(par.b1.m,2,mean)

## MCMC diagnostics

#Calculation of B and W statistics
B.mc1<-rep(0,24)
W.mc1<-rep(0,24)

n=(mcmc-burn)/2
m=chain*2
for (ii in 1:24) {
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

Eff.mc1<-rep(0,24)

for (ii in 1:24) {
	x2<-acf[1:ind[ii],ii]
	Eff.mc1[ii]<-m*n/(1+2*sum(x2))  
	}

Eff.mc1


#Bayesian forecasting based on posterior samples 

nsim<-20000
Mx.LC2.SEA.sim<-matrix(nrow=nsim,ncol=300)
colnames(Mx.LC2.SEA.sim)<-paste(rep(1:30, each=10),AgeGroup)
k.LC2.SEA.sim<-matrix(nrow=nsim,ncol=30)
nr<-nrow(xx1.SEA.m)
kk1.SEA.m<-do.call(rbind,kk1.SEA)
start.x<-colSums(kk1.SEA.m)[length(colSums(kk1.SEA.m))]/nr
k.sim<-c(start.x,rep(100,30))
Mx.sim<-matrix(nrow=10,ncol=30)

set.seed(105)
for(j in 1:nsim){
samp.ind<-sample(x=nr,1,rep=TRUE)
	for (k in 1:30){
	k.sim[k+1]<-k.sim[k]+rnorm(1,mean=xx1.SEA.m[samp.ind,23],
		sd=sqrt(exp(xx1.SEA.m[samp.ind,21])))
	}
k.LC2.SEA.sim[j,]<-k.sim[-1]

mu=xx1.SEA.m[samp.ind,1:10]+xx1.SEA.m[samp.ind,11:20]%*%t(k.sim[-1])

	for (m in 1:30){
	Mx.sim[,m]<-mu[,m]+rnorm(n=10,mean=0,sd=sqrt(exp(xx1.SEA.m[samp.ind,22])))
	}
Mx.LC2.SEA.sim[j,]<-as.vector(exp(Mx.sim))
}

CI95.k.LC2.SEA.sim<-apply(k.LC2.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.k.LC2.SEA.sim<-apply(k.LC2.SEA.sim,2,quantile, probs = c(.005,.5,.995))

CI95.Mx.LC2.SEA.sim<-apply(Mx.LC2.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.Mx.LC2.SEA.sim<-apply(Mx.LC2.SEA.sim,2,quantile, probs = c(.005,.5,.995))






#################### DLM LC with switching k(t) variance 

buildFunB2 <- function(x, vecX,start.x, start.V) {
FF <- matrix(c(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],rep(0,20),x[11],x[12],x[13],
	x[14],x[15],x[16],x[17],x[18],x[19],x[20]),ncol = 4)
GG <- matrix(c(1,0,0,0,0,1,0,1,0,0,1,100,0,0,0,1),ncol=4)
JGG<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0), ncol=4)
W <- diag(c(0,0,0,100))
JW<-diag(c(0,0,0,2))
X<-cbind(rep(0,per.SEA),exp(x[22])*matrix(vecX,ncol=1)+exp(x[21])*(1-matrix(vecX,ncol=1)))
V<-diag(rep(exp(x[23]),10)) 
dlmSEA<-dlm(FF=FF,V=V,GG=GG,JGG=JGG,W=W,JW=JW,X=X,m0=c(1,x[24],x[25],start.x),
	C0=diag(c(0.000001,0.000001,0.000001,start.V))) 
return(dlmSEA)
}

chain<-5
mcmc<-5000
burn<-1000
every<-1

x<-c(rep(-4,10),rep(0.1,14),100)   #Par 25 not used in SEA case
x0<-matrix(nrow=chain,ncol=25)
x0[1,]<-c(rep(-4,10),rep(0.1,14),100)		 
x0[2,]<-c(rep(-4,10),rep(0.1,10),0.05,0.1,0.1,0.1,100)
x0[3,]<-c(rep(-4,10),rep(0.1,10),0.5,1,0.1,0.1,100)
x0[4,]<-c(rep(-4,10),rep(0.1,10),0.1,0.1,1,0.1,100)
x0[5,]<-c(rep(-4,10),rep(0.1,10),0.1,0.1,0.1,1,100)
start.x<-9
start.V<-10

xx2.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=27),simplify=FALSE)
kk2.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=(per.SEA+1)),simplify=FALSE)
Im.SEA<-replicate(chain, matrix(nrow=(mcmc-burn)/every, ncol=per.SEA),simplify=FALSE)
tbl.SEA<-replicate((mcmc-burn)*chain/every, matrix(nrow=2, ncol=2),simplify=FALSE)

p<-matrix(nrow=2, ncol=(per.SEA+1)) 
ind<-matrix(nrow=2, ncol=per.SEA)

prior.beta<-1
prior.beta1<-1  
prior.beta2<-1 
#prior.m<-c(0,0)
prior.m<-0
#prior.m.prec<-matrix(c(1/5,0,0,1/5),ncol=2)
prior.m.prec<-1/5
prior.g1<-2.1
prior.g2<-0.1   
prior.b<-c(-4,0.1) #equal for each age group 
prior.b.prec<-matrix(c(1/5,0,0,1/5),ncol=2)


set.seed(102)
for (ii in 1:chain) {
	x<-x0[ii,]
	i.save<-0
	i0<-0.5  #transition probability from 0 to 0
	i1<-0.5  #transition probability from 1 to 1
	vecX<-rep(0, per.SEA) 

   for (i in 1:mcmc) { 
	
	# 1. Generate sample of states
	mod=buildFunB2(x, vecX, start.x, start.V)
	modFilt <- dlmFilter(t(as.matrix(log(MxM.SEA))), mod, simplify=TRUE)
	k <- dlmBSample(modFilt)  

	
	# 2. Generate mixture index
	# 2.1. generate probabilities
	p[1,1]<-(1-i1)/(2-i1-i0)  
	p[2,1]<-(1-i0)/(2-i1-i0)
	
	for (j in 1:per.SEA) { 
		ind[1,j]<-p[1,j]*i0+p[2,j]*(1-i1)  #prior probabilities of states
		ind[2,j]<-p[2,j]*i1+p[1,j]*(1-i0)  
		
		lk1<-dnorm(x=(k[j+1,4]-k[j,4]),mean=k[j,2],sd=sqrt(exp(x[21]))) 
		lk2<-dnorm(x=(k[j+1,4]-k[j,4]),mean=k[j,2],sd=sqrt(exp(x[22])))

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

	aa<-sqrt(exp(x[22])*matrix(vecX,ncol=1)+exp(x[21])*(1-matrix(vecX,ncol=1))) #vector of sd's 
	xs<-1/aa
	dk<-k[-1,4]-k[-nrow(k),4]
	dks<-dk/aa

	mu0<-rnorm(n=1,mean=(prior.m.prec*prior.m+sum(xs*dks))/(prior.m.prec+sum(xs*xs)), 
		sd=sqrt(1/(prior.m.prec+sum(xs*xs))))


	# 3.3. sample W

	hh=exp(x[22])/exp(x[21])-1    
	dkss<-dk/sqrt(1+hh*vecX)
	ab<-cbind(rep(1,59),c(rep(0,35),rep(1,24)))
	res1<-dkss-mu0/sqrt(1+hh*vecX)
	
	w0.0<-1/rgamma(n=1, prior.g1+per.SEA/2, prior.g2+0.5*sum(res1^2))

	dkss1<-(dk/sqrt(w0.0)-mu0/sqrt(w0.0))*vecX
	res2<-dkss1[dkss1!=0]
	
	hh<-1/rgamma(n=1, prior.g1+length(res2)/2, prior.g2+0.5*sum(res2^2))-1
	
	ifelse(hh>0,
		w1.0<-(1+hh)*w0.0,
		w1.0<-exp(x[22])*w0.0/exp(x[21]))


	# 3.4. sample alphas and betas

	kk<-rbind(rep(1,per.SEA),k[-1,4])%*%cbind(rep(1,per.SEA),k[-1,4])
	b0<-cbind(1:10,1:10)
	for(j in 1:10){	
		
		Mu=solve(prior.b.prec+(1/exp(x[23]))*kk)%*%(prior.b.prec%*%prior.b+(1/exp(x[23]))*
			rbind(rep(1,59),k[-1,4])%*%log(MxM.SEA[j,]))
		Sigma=solve(prior.b.prec+(1/exp(x[23]))*kk)

		b01<-rnorm(n=1,mean=Mu[2],sd=sqrt(Sigma[2,2]))
		b00<-rnorm(n=1,mean=Mu[1]+Sigma[1,2]*(b01-Mu[2])/Sigma[2,2],
			sd=sqrt(Sigma[1,1]-Sigma[1,2]*Sigma[2,1]/Sigma[2,2]))
		b0[j,]<-c(b00,b01)
  	}


	# 3.5. sample V

	res3 <- log(MxM.SEA)-b0[,1]-b0[,2]%*%t(k[-1,4])
	v<-1/rgamma(n=1, prior.g1+per.SEA*10/2, prior.g2+0.5*sum(res3^2))    


	# 3.6. reweight betas and k(t) parameters

	bb<-sum(b0[,2])
	b<-b0[,2]/bb
	a<-b0[,1]+mean(k[-1,4])*b0[,2]
	mu<-(mu0-mean(k[-1,4]))*bb
	w.0<-w0.0*(bb^2)
	w.1<-w1.0*(bb^2)

	# 4. Populate updated paramters to x

	x<-c(a,b,log(w.0),log(w.1),log(v),mu,100)    

	# 5. Save selected itterations

	if(i>burn && !(i %% every)){
		i.save<-i.save+1
		kk2.SEA[[ii]][i.save,]<-k[,4]
		xx2.SEA[[ii]][i.save,]<-c(x,i0,i1)
		Im.SEA[[ii]][i.save,]<-vecX
		tbl.SEA[[(ii-1)*(mcmc-burn)+i.save]]<-tr
	}
   }
}



## Parameter summary

par.b2.m<-matrix(nrow=chain*2, ncol=27)
par.b2.var<-matrix(nrow=chain*2, ncol=27)


for (jj in 1:chain) {
	xx<-xx2.SEA[[jj]]
	xx[,21:23]<-exp(xx[,21:23])
	par.b2.m[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,mean)
	par.b2.m[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,mean)
	par.b2.var[jj*2-1,]<-apply(xx[1:(mcmc-burn)/2,],2,var)
	par.b2.var[jj*2,]<-apply(xx[((mcmc-burn)/2+1):(mcmc-burn),],2,var)
}
#check for significant deviations between results of different half-chains

par.b2<-apply(par.b2.m,2,mean)


## MCMC diagnostics

#Calculation of B and W statistics
B.mc2<-rep(0,27)
W.mc2<-rep(0,27)

n=(mcmc-burn)/2
m=chain*2
for (ii in 1:27) {
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

Eff.mc2<-rep(0,27)

for (ii in 1:27) {
	x2<-acf[1:ind[ii],ii]
	Eff.mc2[ii]<-m*n/(1+2*sum(x2))  
	}

Eff.mc2



## Bayesian forecasting based on posterior samples 

kk2.SEA.m<-do.call(rbind,kk2.SEA)
nsim<-20000
Mx.LC3.SEA.sim<-matrix(nrow=nsim,ncol=300)
colnames(Mx.LC3.SEA.sim)<-paste(rep(1:30, each=10),AgeGroup)
p.sim<-matrix(nrow=2,ncol=31)
vecX.sim<-rep(10,30)
Im.SEA.sim<-matrix(nrow=nsim,ncol=30)
k.LC3.SEA.sim<-matrix(nrow=nsim,ncol=30)
nr<-nrow(xx2.SEA.m)
start.x<-colSums(kk2.SEA.m)[length(colSums(kk2.SEA.m))]/nr
k.sim<-c(start.x,rep(100,30))
Mx.sim<-matrix(nrow=10,ncol=30)


set.seed(107)
for(j in 1:nsim){
	samp.ind<-sample(x=nr,1,rep=TRUE)
	
	#simulate the regime vector
	i0<-xx2.SEA.m[samp.ind,26]
	i1<-xx2.SEA.m[samp.ind,27]
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
	Im.SEA.sim[j,]<-vecX.sim
	
	#simulate k(t)
	var.vect<-vecX.sim*xx2.SEA.m[samp.ind,22]+(1-vecX.sim)*xx2.SEA.m[samp.ind,21]
	for (k in 1:30){
		k.sim[k+1]<-k.sim[k]+rnorm(1,mean=xx2.SEA.m[samp.ind,24],
		sd=sqrt(exp(var.vect[k])))
		}
	k.LC3.SEA.sim[j,]<-k.sim[-1]

	#simulate mortality rates
	mu=xx2.SEA.m[samp.ind,1:10]+xx2.SEA.m[samp.ind,11:20]%*%t(k.sim[-1])
	for (m in 1:30){
		Mx.sim[,m]<-mu[,m]+rnorm(n=10,mean=0,sd=sqrt(exp(xx2.SEA.m[samp.ind,23])))
		}
	Mx.LC3.SEA.sim[j,]<-as.vector(exp(Mx.sim))
	}

CI95.k.LC3.SEA.sim<-apply(k.LC3.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.k.LC3.SEA.sim<-apply(k.LC3.SEA.sim,2,quantile, probs = c(.005,.5,.995))

CI95.Mx.LC3.SEA.sim<-apply(Mx.LC3.SEA.sim,2,quantile, probs = c(.025,.5,.975))
CI99.Mx.LC3.SEA.sim<-apply(Mx.LC3.SEA.sim,2,quantile, probs = c(.005,.5,.995))


## Save data for VAR calculations

write.table(CI99.Mx.LC3.SEA.sim,".\\Results\\CI99_Mx_LC3_SEA_sim.txt",sep="\t")
write.table(Mx.LC3.SEA.sim,".\\Results\\Mx_LC3_SEA_sim.txt",sep="\t", row.names=FALSE)
write.table(par.b2,".\\Results\\par_b2.txt",sep="\t",row.names=FALSE)
write.table(start.x,".\\Results\\start_x.txt",sep="\t",row.names=FALSE)
write.table(k.LC3.SEA.sim,".\\Results\\k_LC3_SEA_sim.txt",sep="\t", row.names=FALSE)

## Save data for conditional likelihood simulations

write.table(par.b2,".\\Results\\par_b2.txt",sep="\t",row.names=FALSE)
write.table(kk2.SEA.m,".\\Results\\kk2_SEA_m.txt",sep="\t",row.names=FALSE)
write.table(MxM.SEA,".\\Results\\MxM_SEA.txt",sep="\t", row.names=FALSE)


