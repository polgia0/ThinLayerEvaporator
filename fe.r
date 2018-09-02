# activity coefficient :UNIQUAC
uniquac<-function(x,T){
	z<-10
	rr<-c(2.57,0.92)
	qq<-c(2.34,1.40)
	aa<-matrix(c(0,-100.71,530.99,0),2,2)
	taum<-matrix(rep(0,4),2,2)
	fi<-x*rr/sum(x*rr)
	theta<-x*qq/sum(x*qq)
	taum<-exp(-aa/T)
	sumtau<-apply(theta*taum,2,sum)
	sum2<-apply(theta*t(taum)/sumtau,2,sum)
	L<-z/2*(rr-qq)-(rr-1)
	sum1x<-sum(L*x)
	lngam<-log(fi/x)+z/2*qq*log(theta/fi)+L-fi/x*sum1x-qq*log(sumtau)+qq-qq*sum2
    gam<-exp(lngam)
return(gam)}
# tensione vapore acqua in bar temperatura in k
P0W<-function(T){
	P0W<-10^(29.8605-3.1522e3/T-7.3037*log10(T)+2.4247e-9*T+1.8090e-6*T^2)  #mmHg
	P0W<-P0W/760                                  #bar
return(P0W)}
# tensione vapore acetone in bar temperatura in k
P0Ac<-function(T){
	P0Ac<-10^(28.5884-2.4690e3/T-7.3510*log10(T)+2.8025e-10*T+2.7361e-6*T^2)#mmHg
	P0Ac<-P0Ac/760                                #bar
return(P0Ac)}
do.transport<-function(x1,pars){
    with (as.list(pars),{
        x2<-1-x1
        w1<-MW1*x1/(MW1*x1+MW2*x2)
        w2<-1-w1
        rhomass<-1/(w1/rho1+w2/rho2)          #kg/m3
        Cp<-Cp1*w1+Cp2*w2                     #kcal/kgK
        tc<-1/nrot/Bn                         #s
        h<-2*sqrt(rhomass*Kc*Cp/pi/tc)        #kcal/smK
        rhomol<-1/(x1*MW1/rho1+x2*MW2/rho2)   #kmol/m3
        kc<-2*sqrt(D12/mu/pi/tc)              #m/s
        kx<-kc*rhomol                         #kmol/m2s
        return(list(h=h,kx=kx))
    })
}
reseb<-function(Teb,pars){
  with (as.list(c(pars)),{
    x2<-1-x1
    gam<-uniquac(c(x1,x2),Teb)
    reseb<-P0Ac(Teb)*gam[1]*x1/P+P0W(Teb)*gam[2]*x2/P-1
    return(reseb)
  })  
}
algebric_system<-function(x,pars){
  with (as.list(c(x,pars)),{
	tr<-do.transport(x1,pars)
	kx<-tr$kx
	h<-tr$h/100
	#print(c(kx,h))
	x1s<-x1
	x2s<-1-x1
	for (i in 1:100){
		xsold<-x1s
		# T1<-reseb(223,pars)
		# T3<-reseb(473,pars)
		# print(c(T1,T3))
		Tl<-uniroot(reseb,interval=c(223,473),c(x,pars))
		T2<-Tl$root
		# T1<-223
		# T3<-473
		# gam<-uniquac(c(x1s,x2s),T1)
		# F1<-P0Ac(T1)*gam[1]*x1s/P+P0W(T1)*gam[2]*x2s/P-1
		# gam<-uniquac(c(x1s,x2s),T3)
		# F3<-P0Ac(T3)*gam[1]*x1s/P+P0W(T3)*gam[2]*x2s/P-1
		# for (j in 1:100){
		# 	T2<-(T1+T3)/2
		# 	err<-abs(T2-T1)/T2*100
		# 	if(err<=1e-4)break
		# 	gam<-uniquac(c(x1s,x2s),T2)
		# 	F2<-P0Ac(T2)*gam[1]*x1s/P+P0W(T2)*gam[2]*x2s/P-1
		# 	FF<-F1*F2
		# 	if(FF>0){
		# 		T1<-T2
		# 		F1<-F2
		# 	}else{
		# 		T3<-T2
		# 		F3<-F2
		# 	}
		# }
		x2s<-1-x1s
		gam<-uniquac(c(x1s,x2s),T2)
		y1s<-P0Ac(T2)*gam[1]*x1s/P
		y2s<-P0W(T2)*gam[2]*x2s/P
		Nt<-h*(TH-T2)/((DH1-DH2)*y1s+DH2)                   #kmol/sm
		x1s<-(Nt/kx+1)*x1-Nt*y1s/kx                         #m
		err<-abs((x1s-xsold)/x1s)*100
		if(err<0.1) break
		#print(c(x1s,x2s,y1s,y2s,T2))
	}
	Q<-h*(TH-T2)*C
	N1<-Nt*y1s
	return(c(x1eq=x1s,x2eq=x2s,y1eq=y1s,y2eq=y2s,T=T2,Q=Q,Nt=Nt,N1=N1))
	})
}
Differential_system<-function(t,y,pars){
with(as.list(c(y,pars)),{
    x<-c(x1eq=x1,x2eq=1-x1,y1eq=1,y2eq=0,T=TH,Q=0,Nt=0,N1=0)
    root<-algebric_system(x,c(y,pars))
	  names(root)<-names(x)
    with(as.list(root),{
		ydot<-rep(0,length(y))
		ydot[1]<--Nt*C
		ydot[2]<-Nt*C/L*(x1-y1eq)
		ydot[3]<-Q
		ydot[4]<-Nt*C
		ydot[5]<-N1*C
		return(list(ydot,root))
    })
})
}
Constrains<-function(t,y,pars)
{
	yroot<-rep(1,2)
    with(as.list(c(y,pars)),{
		yroot[1]<-L-0.0
		yroot[2]<-x1-0.0
        return(yroot)
    })
}
eventfun<-function(t,y,pars)
{
    yroot<-Constrains(t,y,pars)
	indroot<-as.logical(abs(yroot)<1e-9)
    with(as.list(c(y,pars)),{
	    if(indroot[1]){
		    print('Any Liquid Flow')
		}
		if(indroot[2]){
		    print('Any Acetone Flow')
		}
		return(y)
	})
}
require(gdata)
require(deSolve)
require(rootSolve)
pars<-c(
	MW1=58.00,                                                      #kg/kmol
	MW2=18.00,                                                      #kg/kmol
	rho1=770,                                                         #kg/m3
	rho2=970,                                                         #kg/m3
	Cp1=0.517,                                                     #kcal/kgK
	Cp2=1.004,                                                     #kcal/kgK
	Dm1=1.14e-9,                                                       #m2/s
	Dm2=2.55e-9,                                                       #m2/s
	D12=1.50,                                                          #m2/s
	DH1=7481,                                                     #kcal/kmol
	DH2=9799,                                                     #kcal/kmol
	mu=100,                                                               #-
	Kc=1.1e-4,                                                     #kcal/smK
	W0=0.375,                                                             #-
	LFw=350,                                                           #kg/h
	WF=0.15,                                                              #-
	TH=55,                                                                #C
	P=0.15,                                                             #bar
	nrot=425,                                                           #rpm
	Bn=3,                                                                 #-
	Z=0.73,                                                               #m
	D=0.22)                                                               #m
pars<-unlist(
    within(as.list(pars),{
		LFw<-LFw/3600  
		nrot<-nrot/60	
		TH<-TH+273
	})
)
pars<-with(as.list(pars),
	c(pars,
	C=pi*D,
    xF=WF/MW1/(WF/MW1+(1-W0-WF)/MW2),
    LF=LFw*(WF/MW1+(1-W0-WF)/MW2)
    )
)
attach(as.list(pars))
dz<-log(1:60)
z<-rev(pars["Z"]*(dz[60]-dz)/dz[60])
y<-with(c(as.list(pars)),{                                 # initial conditions
c(  L=LF,                                                                #mol/s
    x1=xF,                                                                   #-
	Qt=0,                                                               #kcal/s
	DIS=0,                                                                  #kg/s
	DIS1=0)                                                                    #-
})
out<-lsodar(
	y=y,
	times=z,
	fun=Differential_system,
	rootfun=Constrains,
	rtol=1e-4,
	atol=rep(1e-6,length(y)),
	parms=pars,
	events=list(func=eventfun,root=TRUE,terminalroot=c(1,2))
	)
out<-as.data.frame(out)
out<-with(as.list(pars),{
within(out,{
    T<-T-273.15                                              #conversion K to C
})
})
  with(c(as.list(pars)),{ 
  	print('----------------------------------',quote=F)
  	print(paste('Feed (kg/h)    =',format(LFw,digits=4),sep=''),quote=F)
  	print(paste('WF acetone     =',format(WF*100,digits=4),sep=''),quote=F)
  	print(paste('Rotation (rpm) =',format(nrot*60,digits=4),sep=''),quote=F)
  	print(paste('Temperature (C)=',format(TH,digits=4),sep=''),quote=F)
  	print(paste('Pressure (bar) =',format(P,digits=4),sep=''),quote=F)
  	print('----------------------------------',quote=F)
  	nr<-nrow(out)
  	print(out)
  Tcond<-out$T[nr]
  MWave<-MW1*out$x1[nr]+MW2*(1-out$x1[nr])
  Lwout<-out$L[nr]*MWave+LFw*W0
  ER<-(LFw-Lwout)/LFw*100
  ERpb<-(LF-out$L[nr])/LF*100
  Ex1<-(LF*xF-out$L[nr]*out$x1[nr])/LF/xF*100
  Ex2<-(LF*(1-xF)-out$L[nr]*(1-out$x1[nr]))/LF/(1-xF)*100
  Pw<-out$Qt[nr]/C/Z*4.183*10
  xfin<-out$x1[nr]
  Lfin<-out$L[nr]*(xfin*MW1+(1-xfin)*MW2)+LFw*W0
  WFin<-out$L[nr]*xfin*MW1/Lfin*100
  print(paste('Total Evaporation(%) =',format(ER,digits=4),sep=''),quote=F)
  print(paste('Acetone (%)          =',format(Ex1,digit=4),sep=''),quote=F)
  print(paste('Water   (%)          =',format(Ex2,digit=4),sep=''),quote=F)
  print(paste('W end   (-)          =',format(WFin,digit=4),sep=''),quote=F)
print(paste('Thermal Power(kW/m2) =',format(Pw,digit=4),sep=''),quote=F)
print(paste('T cond (C)           =',format(Tcond,digit=4),sep=''),quote=F)
print('----------------------------------',quote=F)
})
  