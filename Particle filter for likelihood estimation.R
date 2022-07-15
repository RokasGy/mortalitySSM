##### Open stochastic mortaly projections files

#Set working directory for input/output
#SPECIFY YOUR WORKING DIRECTORY HERE
#setwd("C:\\....\\")

per.SEA<-118

# Estimated parameters
par.b2<- read.table(".\\Results\\par_b2.txt", header=TRUE, row.names=NULL)
par.b2<-t(par.b2)

# Fitted historic k(t) parameters
kk2.SEA.m<- read.table(".\\Results\\kk2_SEA_m.txt", 
check.names=FALSE, header=FALSE, row.names=NULL)
kk2.SEA.m<-as.matrix(kk2.SEA.m)

# Historic mortality rates
MxM.SEA<- read.table(".\\Results\\MxM_SEA.txt", 
check.names=FALSE, header=TRUE, row.names=NULL)
MxM.SEA<-as.matrix(MxM.SEA)


##### Conditional likelihood simulation

mcpf<-5000   #the the number of "particles"
i0.pf<-par.b2[26]  
i1.pf<-par.b2[27]
m.pf<-par.b2[24] # one drift parameter for SE
w0.pf<-par.b2[21]
w1.pf<-par.b2[22]
v.pf<-par.b2[23]
a.pf<-par.b2[1:10]   
b.pf<-par.b2[11:20]
reg.pf<-rep(list(matrix(nrow=mcpf, ncol=3)),per.SEA+1)
# separate list for each year of fitted data

lk.pf<-matrix(nrow=mcpf,ncol=per.SEA)

set.seed(108)
#initiallize particle filter
reg.pf[[1]]<-cbind(sample(x=kk2.SEA.m[,1], size=mcpf, replace=TRUE),rbinom(n=mcpf,prob=(1-i0.pf)
	/(2-i0.pf-i1.pf),size=1),rep(10,mcpf)) 

for (j in 1:per.SEA){	
	#one step ahead forecasts
	xm.pf<-reg.pf[[j]][,1]+m.pf	#one step ahead forecast of k(t)
	xw0.pf<-reg.pf[[j]][,3]+w0.pf  #if S(t)=0
	xw1.pf<-reg.pf[[j]][,3]+w1.pf  #if S(t)=1

	#residuals and expected V variance
	res.pf<-apply(as.matrix(xm.pf,ncol=1,nrow=mcpf),1,function(x) log(MxM.SEA[,j])-a.pf-b.pf*x)  ##3
	vx0.pf<-lapply(xw0.pf,function(x) as.matrix(b.pf,ncol=1)%*%x%*%t(b.pf)+diag(v.pf,10))
	vx0.pf.inv<-lapply(vx0.pf, function (x) solve(x))
	vx1.pf<-lapply(xw1.pf,function(x) as.matrix(b.pf,ncol=1)%*%x%*%t(b.pf)+diag(v.pf,10))
	vx1.pf.inv<-lapply(vx1.pf, function (x) solve(x))

	#calculate weights of regimes 0 and 1
	trActto0<-reg.pf[[j]][,2]*(1-i1.pf)+(1-reg.pf[[j]][,2])*i0.pf
	trActto1<-reg.pf[[j]][,2]*i1.pf+(1-reg.pf[[j]][,2])*(1-i0.pf)
	s.det0<-sapply(vx0.pf, function(x) 1/sqrt(det(x)))
	s.det1<-sapply(vx1.pf, function(x) 1/sqrt(det(x)))
	eta0.pf<-sapply(1:mcpf,function(i) s.det0[i]*exp(-0.5*(t(res.pf[,i])
		%*%vx0.pf.inv[[i]]%*%res.pf[,i])))
	eta1.pf<-sapply(1:mcpf,function(i) s.det1[i]*exp(-0.5*(t(res.pf[,i])
		%*%vx1.pf.inv[[i]]%*%res.pf[,i])))
	eta0.pf<-eta0.pf*trActto0
	eta1.pf<-eta1.pf*trActto1
	
	#sample S(t|y)
	pr0.pf<-sapply(1:mcpf,function(i) eta0.pf[i]/(eta0.pf[i]+eta1.pf[i]))
	smpS.pf<-sapply(pr0.pf, function(x) rbinom(n=1, size=1,prob=1-x))	

	#sample x(t|y,S)
	xw.pf<-(1-smpS.pf)*xw0.pf+smpS.pf*xw1.pf
	vx.pf<-lapply(xw.pf,function(x) as.matrix(b.pf,ncol=1)%*%x%*%t(b.pf)+diag(v.pf,10))
	vx.pf.inv<-lapply(vx.pf, function (x) solve(x))
	det<-sapply(vx.pf, function(x) 1/sqrt(det(x)))

	ym.pf<-sapply(1:mcpf,function(i) xm.pf[i]+xw.pf[i]%*%t(as.matrix(b.pf))%*%
		vx.pf.inv[[i]]%*%res.pf[,i])
	yw.pf<-sapply(1:mcpf,function(i) xw.pf[i]-xw.pf[i]%*%t(as.matrix(b.pf))%*%
		vx.pf.inv[[i]]%*%as.matrix(b.pf)%*%xw.pf[i])
	smpX.pf<-sapply(1:mcpf, function (i) rnorm(n=1, mean=ym.pf[i], sd=sqrt(yw.pf[i])))

	#update vector of particles
	reg.pf[[j+1]][,1]<-smpX.pf
	reg.pf[[j+1]][,2]<-smpS.pf
	reg.pf[[j+1]][,3]<-yw.pf

	#likelihood calculation
		
	lk.pf[,j]<-sapply(1:mcpf,function(i) {
	trActto0[i]*s.det0[i]*exp(-0.5*(t(res.pf[,i])%*%vx0.pf.inv[[i]]%*%res.pf[,i]))+
	trActto1[i]*s.det1[i]*exp(-0.5*(t(res.pf[,i])%*%vx1.pf.inv[[i]]%*%res.pf[,i]))
		})
	}	

# Results:
per.SEA*log((2*pi)^(-5))+sum(log(colSums(lk.pf)/mcpf))


