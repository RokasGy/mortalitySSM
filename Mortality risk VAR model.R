##### Open stochastic mortaly projections files

#Set working directory for input/output
#SPECIFY YOUR WORKING DIRECTORY HERE
#setwd("C:\\....\\")

# 99% confidence intervals
CI99.Mx.LC3.SEA.sim<- read.table(".\\Results\\CI99_Mx_LC3_SEA_sim.txt", 
check.names=FALSE)

# stochastic mortality simulations
Mx.LC3.SEA.sim<- read.table(".\\Results\\Mx_LC3_SEA_sim.txt", 
check.names=FALSE, header=TRUE, row.names=NULL)
Mx.LC3.SEA.sim<-as.matrix(Mx.LC3.SEA.sim)

# Estimated parameters
par.b2<- read.table(".\\Results\\par_b2.txt", header=TRUE, row.names=NULL)

# Starting point for k(t) projections
start.x<- read.table(".\\Results\\start_x.txt", header=TRUE, row.names=NULL)
start.x<-start.x[1,1]

# k(t) simulated projections
k.LC3.SEA.sim<- read.table(".\\Results\\k_LC3_SEA_sim.txt", 
check.names=FALSE, header=TRUE, row.names=NULL)
k.LC3.SEA.sim<-as.matrix(k.LC3.SEA.sim)




##### Run-off VaR calculation

# base references vector and output matrixes
colref.base<-c(seq(from=1, by=10, length.out=5),seq(from=52, by=10, length.out=5),
	seq(from=103, by=10, length.out=5),seq(from=154, by=10, length.out=5),
	seq(from=205, by=10, length.out=5),seq(from=256, by=10, length.out=5))

col.names<-c("2530", "3030", "3530","4030","4530","2520", "3020","3520","4020","4520","5020","5520",
	"2510", "3010", "3510","4010","4510","5010","5510","6010","6510")
q.SEA.lev<-matrix(nrow=20000,ncol=(9+7+5))
colnames(q.SEA.lev) <- col.names
q.SEA.lev.CI<-matrix(nrow=2,ncol=(9+7+5))
colnames(q.SEA.lev.CI) <- col.names
VaR1.SEA.lev<-matrix(nrow=1,ncol=(9+7+5))
colnames(VaR1.SEA.lev) <- col.names

q.SEA.dec<-matrix(nrow=20000,ncol=(9+7+5))
colnames(q.SEA.dec) <- col.names
q.SEA.dec.CI<-matrix(nrow=2,ncol=(9+7+5))
colnames(q.SEA.dec.CI) <- col.names
VaR1.SEA.dec<-matrix(nrow=1,ncol=(9+7+5))
colnames(VaR1.SEA.dec) <- col.names

Qx.LC3.SEA.f<-matrix(1-exp(-CI99.Mx.LC3.SEA.sim[2,]),nrow=1)
Px.LC3.SEA.sim<-exp(-Mx.LC3.SEA.sim)
Qx.LC3.SEA.f<-as.numeric(Qx.LC3.SEA.f)

for (t in c(30, 20, 10)) {

	max.age=11-t/5
	for (a in 1:max.age){ 
	
# Level benefits
		aa<-1-apply(Px.LC3.SEA.sim[,head(colref.base+(a-1),t)],1,prod)
		q.SEA.lev[,paste0(25+(a-1)*5, t)]<-aa
		aa.CI<-quantile(aa, probs = c(.5,.995))
		q.SEA.lev.CI[,paste0(25+(a-1)*5, t)]<-aa.CI
	
		f.sf<-function(sf) {	#function for calculation of VaR rate (level benefits)
			aa.sf<-1-Qx.LC3.SEA.f*(1+sf)
			diff.sf<-abs(aa.CI[2]-(1-prod(aa.sf[head(colref.base+(a-1),t)])))
			}
		
		aa.rate<-optimize(f.sf, maximum=FALSE, c(-3, 3))
		#the function searches for the level shift in mortality that matches sochastically 
		#calcululated confidence level

		VaR1.SEA.lev[,paste0(25+(a-1)*5, t)]<-aa.rate$minimum

# Decreasing benefits
		bb.list<-list()
			for (j in 1:t){
				bb.list[[j]]<-head(colref.base+(a-1),t)[1:(t+1-j)]
			}  	#create references for calculation of mpx for different m
	
		q.sim<-matrix(nrow=20000,ncol=t) 
			for (i in 1:(t-1)){
				q.sim[,(t+1-i)]<-1-apply(Px.LC3.SEA.sim[,bb.list[[i]]],1,prod)
			}
			q.sim[,1]<-1-Px.LC3.SEA.sim[,bb.list[[t]]] #first row - no need for multiplications
	
		bb<-(1/t)*apply(q.sim,1,sum)
		q.SEA.dec[,paste0(25+(a-1)*5, t)]<-bb
		bb.CI<-quantile(bb,probs = c(.5,.995))
		q.SEA.dec.CI[,paste0(25+(a-1)*5, t)]<-bb.CI 
		
		f.sf.d<-function(sf) {
			bb.sf<-1-Qx.LC3.SEA.f*(1+sf)
			qbb.sf<-matrix(nrow=1,ncol=t)
			for (i in 1:(t-1)){
				qbb.sf[(t+1-i)]<-1-prod(bb.sf[bb.list[[i]]])
			}
			qbb.sf[,1]<-1-bb.sf[bb.list[[t]]]
			sum.sf<-(1/t)*sum(qbb.sf)
			diff.sf<-abs(bb.CI[2]-sum.sf)
		}
	
		bb.rate<-optimize(f.sf.d, maximum=FALSE, c(-2, 2))
		VaR1.SEA.dec[,paste0(25+(a-1)*5, t)]<-bb.rate$minimum
	}
}


# Results:
# the first of two digit numbers (in column heading) indicates age at policy inception, the second - policy term
# VAR rates are expressed as one-off level shock on base estimate mortality rates required
# to achieve the 99.5% confidence interval of the Best Estimate Liability. 

# Level benefits
VaR1.SEA.lev

# Decreasing benefits
VaR1.SEA.dec



## 1 year VaR

#base references vector and output matrixes
colref.base2<-c(seq(from=1, by=10, length.out=4),seq(from=42, by=10, length.out=5),
seq(from=93, by=10, length.out=5),seq(from=144, by=10, length.out=5),
seq(from=195, by=10, length.out=5),seq(from=246, by=10, length.out=5))

col.names<-c("2530", "3030", "3530","4030","4530","2520", "3020","3520","4020","4520","5020","5520",
"2510", "3010", "3510","4010","4510","5010","5510","6010","6510")
kka<-matrix(nrow=20000,ncol=(9+7+5))
colnames(kka) <- col.names
kkb<-matrix(nrow=2,ncol=(9+7+5))
colnames(kkb) <- col.names
kkd<-matrix(nrow=1,ncol=(9+7+5))
colnames(kkd) <- col.names

kke<-matrix(nrow=20000,ncol=(9+7+5))
colnames(kke) <- col.names
kkf<-matrix(nrow=2,ncol=(9+7+5))
colnames(kkf) <- col.names
kkh<-matrix(nrow=1,ncol=(9+7+5))
colnames(kkh) <- col.names
kki = matrix(nrow=20000, ncol=29)



# VaR calculations
for (d in c(0.05, 0.1)){ #calculations for 2 reserve sensitivities
# adjust sensitivities if needed

	# simulated future trend/ drift parameter
	dd<-par.b2[24]+d*(k.LC3.SEA.sim[,1]-start.x-par.b2[24]) 

	# calculation of projected mean k(t) for t+1, t+2, etc.
	for (i in 1:20000){
		kki[i,]<-k.LC3.SEA.sim[i,1]+seq(1:29)*dd[i]
	}

	# reestimation of periodic yearly Mx based on simulated k(t), assuming a(x) and b(x) are fixed
	kkl<-apply(kki, 1, function(x) exp(par.b2[11:20]%*%t(x)+par.b2[1:10]))  
	kkl<-t(kkl) #transpose the matrix that rows represent simulations
	colnames(kkl)<-paste(rep(2:30, each=10),AgeGroup) 
	kkm<-exp(-kkl) 	#Calculate Px


	for (t in c(30, 20, 10)) { #calculations for 3 different terms to maturity
	
		max.age=11-t/5
		for (a in 1:max.age){ #calculations for 3 age groups starting from 25 years

#level benefits
		cc.BE<-1-apply(kkm[,head(colref.base2+(a-1),t-1)],1,prod) #simulated BE's assuming that survived first year
		cc<-(1-Px.LC3.SEA.sim[,1+(a-1)])+Px.LC3.SEA.sim[,1+(a-1)]*cc.BE   #first year loss plus BE
		kka[,paste0(25+(a-1)*5, t)]<-cc
		cc.CI<-quantile(cc,probs = c(.5,.995))
		kkb[,paste0(25+(a-1)*5, t)]<-cc.CI
		f.sf<-function(sf) {
			cc.sf<-1-Qx.LC3.SEA.f*(1+sf)
			diff.sf<-abs(cc.CI[2]-(1-prod(cc.sf[head(colref.base+(a-1),t)])))
		}
		cc.rate<-optimize(f.sf, maximum=FALSE, c(-2, 2))
		kkd[,paste0(25+(a-1)*5, t)]<-cc.rate$minimum 

#decreasing benefits

		dd.list<-list()
		for (j in 1:(t-1)){   
			dd.list[[j]]<-head(colref.base2+(a-1),t-1)[1:(t-j)]  
		}  
	
		bb.list<-list()
		for (j in 1:t){
			bb.list[[j]]<-head(colref.base+(a-1),t)[1:(t+1-j)]
		}  

		q.sim<-matrix(nrow=20000,ncol=t-1) #matrix of simulated s year survival probabilities
		for (i in 1:(t-2)){
			q.sim[,(t-i)]<-1-apply(kkm[,dd.list[[i]]],1,prod)
		}
		q.sim[,1]<-1-kkm[,dd.list[[t-1]]] 
	
		dd.BE<-(1/t)*apply(q.sim,1,sum)
		dd<-(1-Px.LC3.SEA.sim[,1+(a-1)])+Px.LC3.SEA.sim[,1+(a-1)]*dd.BE
 		kke[,paste0(25+(a-1)*5, t)]<-dd 
		dd.CI<-quantile(dd,probs = c(.5,.995))
		kkf[,paste0(25+(a-1)*5, t)]<-dd.CI

		f.sf.d<-function(sf) {        
			dd.sf<-1-Qx.LC3.SEA.f*(1+sf)
			p.dd.sf<-matrix(nrow=1,ncol=t)
			for (i in 1:(t-1)){
				p.dd.sf[,(t+1-i)]<-prod(dd.sf[bb.list[[i]]])
			} 
			p.dd.sf[,1]<-dd.sf[bb.list[[t]]]
			q.dd.sf<-1-p.dd.sf
			CI.dd.sf<-(1/t)*sum(q.dd.sf)
			diff.sf<-abs(dd.CI[2]-CI.dd.sf)
		}  
	
		dd.rate<-optimize(f.sf.d, maximum=FALSE, c(-2, 2))
		kkh[,paste0(25+(a-1)*5, t)]<-dd.rate$minimum
		}
	}

assign(paste0("VaR1.SEA.1y.lev.",d*100),kkd)

assign(paste0("VaR1.SEA.1y.dec.",d*100),kkh)
}

# Results:
# the first of two digit numbers (in column heading) indicates age at policy inception, the second - policy term
# VAR rates are expressed as one-off level shock on base estimate mortality rates required
# to achieve the 99.5% confidence interval of the Best Estimate Liability. 

# Level benefits
VaR1.SEA.1y.lev.5 # 5% reserve sensitivity
VaR1.SEA.1y.lev.10 # 10% reserve sensitivity

# Decreasing benefits
VaR1.SEA.1y.dec.5
VaR1.SEA.1y.dec.10

