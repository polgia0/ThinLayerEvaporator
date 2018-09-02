rm(list = ls())
P01<-function(T,pars){
  P01<-10^(pars["P011"]+pars["P012"]/T+pars["P013"]*log10(T)+pars["P014"]*T+pars["P015"]*T^2) #mmHg
  P01<-P01/760  #bar
  return(P01)}
P02<-function(T,pars){
  P02<-10^(pars["P021"]+pars["P022"]/T+pars["P023"]*log10(T)+pars["P024"]*T+pars["P025"]*T^2) #mmHg
  P02<-P02/760  #bar
  return(P02)}
CP1<-function(T,pars){
  CP1<-pars["CPA1"]+pars["CPB1"]*T+pars["CPC1"]*T^2+pars["CPD1"]*T^3    #J/molK
  CP1<-CP1/4.18/pars["MW1"]                     #Kcal/kgK
}
CP2<-function(T,pars){
  CP2<-pars["CPA2"]+pars["CPB2"]*T+pars["CPC2"]*T^2+pars["CPD2"]*T^3    #J/molK
  CP2<-CP2/4.18/pars["MW2"]                     #Kcal/kgK
}
RO1<-function(T,pars){
  RO1<-pars["roA1"]*pars["roB1"]^(-(1-T/pars["roTc1"])^pars["ron1"])    #kg/l
  RO1<-RO1*1000                     #kg/m3
}
RO2<-function(T,pars){
  RO2<-pars["roA2"]*pars["roB2"]^(-(1-T/pars["roTc2"])^pars["ron2"])    #kg/l
  RO2<-RO2*1000                     #kg/m3
}
# activity coefficient :UNIQUAC
uniquac<-function(x,T,pars){
  z<-10
  rr<-c(pars["r1"],pars["r2"])
  qq<-c(pars["q1"],pars["q2"])
  aa<-matrix(c(0,pars["a12"],pars["a21"],0),2,2)
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
  names(gam)<-c("g1","g2")
  return(gam)}
reseb<-function(Tx,pars){
  with (as.list(c(pars)),{
    x2eq<-1-x1eq
    if(x2eq!=0){
        gam<-uniquac(c(x1eq,x2eq),Tx,pars)
        reseb<-P01(Tx,pars)*gam[1]*x1eq/P+P02(Tx,pars)*gam[2]*x2eq/P-1
    }else{
        reseb<-P01(Tx,pars)/P-1
    }
    return(reseb)
  })  
}
resxeq<-function(x1eq,pars){
  with (as.list(c(pars)),{
    Tl<-uniroot(reseb,interval=c(223,473),c(x1eq=as.numeric(x1eq),pars))
    Teb<-Tl$root
    x2eq<-1-x1eq
    gam<-uniquac(c(x1eq,x2eq),Teb,pars)
    Keq<-P01(Teb,pars)*gam[1]/P
    y1eq<-Keq*x1eq
    y2eq<-1-y1eq
    Nt<-h*(TH-Teb)/(DH1*y1eq+DH2*y2eq) #kmol/sm2
    return(x1eq-(Nt+kx)*x1/(Nt*Keq+kx))
  })  
}
algebric_system<-function(y,pars){
  with (as.list(c(y,pars)),{
    if(L>0){
        ws<-L*x1*MW1                          #kg/s
        wns<-L*(1-x1)*MW2                     #kg/s
        wt<-ws+wns
        ws<-ws/wt
        wns<-wns/wt
        Kcm<-Kc1*ws
        Rom<-RO1(Teb,pars)*ws
        Cpm<-CP1(Teb,pars)*ws
        if(wns>0){
          Kcm<-Kcm+Kc2*wns
          Rom<-Rom+RO2(Teb,pars)*wns
          Cpm<-Cpm+CP2(Teb,pars)*wns
        }
        tc<-1/nrot/Bn                         #s
        Re1<-4*LFw/Mu
        Re2<-D^2*nrot*Ro/Mu
        Pr<-Mu*Cpm/Kcm
        #ho<-2*sqrt(Rom*Kcm*Cpm/pi/tc)/he       #kcal/sm2K
        #h<-Kcm/D*(0.018*(Re1^0.46)*(Re2^0.6)*(Pr^0.87)*(D/Z)^0.48*nrot^0.24)
        h<-he1*nrot^he2*LFw^he3    
        if(h<0)h<-0
        kx<-2*Rom*sqrt(D1/pi/tc)              #kg/m2s
        Qv<-LFw/Rom                           #m3/s
        ReL<-Qv/pi/D/Mu                       #-
        wf<-0.33*sqrt(ReL)                    #m/s
        wm<-0.0165*(ReL)**(1/3)               #m/s
        Hv<-Qv/wf/pi/D+0.6*H0*(1-wm/wf)       #m
        Qf<-pi*D*(Hv-0.6*H0)*wf
        Hf<-sqrt(4/pi*Qf/Bn/wf)
        Qm<-Qv-Qf
        Hm<-Qm/pi/D/wm
        if(x1!=1){ # case of binary mixture with water
            xl<-uniroot(resxeq,interval=c(1e-15,0.99999),c(x1=as.numeric(x1),
                                                         h=as.numeric(h),
                                                         kx=as.numeric(kx),
                                                         tc=as.numeric(tc),
                                                         pars))
            x1eq<-xl$root
            x2eq<-1-x1eq
            Tl<-uniroot(reseb,interval=c(223,473),c(x1eq=as.numeric(x1eq),pars))
            Teb<-Tl$root
            gam<-uniquac(c(x1eq,x2eq),Teb,pars)
            y1eq<-P01(Teb,pars)*gam[1]*x1eq/P
            y2eq<-P02(Teb,pars)*gam[2]*x2eq/P
        }else{
            x1eq<-1
            x2eq<-0
            Tl<-uniroot(reseb,interval=c(223,473),c(x1eq=as.numeric(x1eq),pars))
            Teb<-Tl$root
            y1eq<-1
            y2eq<-0
        }
        Nt<-h*(TH-Teb)/(DH1*y1eq+DH2*y2eq) #kmol/sm2
        Q<-h*(TH-Teb)*C
        N1<-Nt*y1eq
        #print(c(h,he1,he2,he3,he4,he5,he6))
        return(c(x1eq=x1eq,x2eq=x2eq,y1eq=y1eq,y2eq=y2eq,Teb=Teb,Q=Q,Nt=Nt,N1=N1,
                 h=h,Qv=Qv,Qf=Qf,Qm=Qm,wf=wf,wm=wm,Hv=Hv,Hm=Hm,Hf=Hf,ReL=ReL,kx=kx))
    }else{
      return(c(x1eq=0,x2eq=0,y1eq=0,y2eq=0,Teb=TH,Q=0,Nt=0,N1=0,
               h=0,Qv=0,Qf=0,Qm=0,wf=0,wm=0,Hv=0,Hm=0,Hf=0,ReL=0,kx=0))
    }
  })
}
Differential_system<-function(t,y,pars){
  with(as.list(c(t,y,pars)),{
    x<-c(x1eq=x1,x2eq=1-x1,y1eq=1,y2eq=0,Teb=TH,Q=0,Nt=0,N1=0,h=0,
         Qv=0,Qf=0,Qm=0,wf=0,wm=0,Hv=0,Hm=0,Hf=0,ReL=0,kx=0)
    root<-algebric_system(x,c(z=t,y,pars))
    names(root)<-names(x)
    with(as.list(root),{
      ydot<-rep(0,length(y))
      ydot[1]<--Nt*C
      ydot[2]<-Nt*C/L*(x1-y1eq)
      ydot[3]<-Q
      ydot[4]<-Nt*C*(y1eq*MW1+y2eq*MW2)
      ydot[5]<-N1*C*y1eq*MW1
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
      print('Any Solvent Flow')
    }
    return(y)
  })
}
server = function(input, output) { # begin server
  values<-reactiveValues(df=NULL,dfopt=NULL,DS=NULL,pars=NULL,ini=NULL,Z=NULL,D=NULL,Bn=NULL,H0=NULL,
                         Rot=NULL,P=NULL,TH=NULL,LFw=NULL,he=NULL,Kc=NULL,Ro=NULL,Mu=NULL,W0=NULL,
                         WF=NULL,NSS=NULL,SS=NULL,HEopt=NULL)
  values$df<-read.csv("compound_data.csv",header=TRUE,sep =";",row.names=1,dec=",",stringsAsFactors=FALSE)
  values$ini<-read.csv("examples_data.csv",header=TRUE,sep =";",row.names=1,dec=",",stringsAsFactors=FALSE)
  output$uiZ<-renderUI({
    numericInput("nZ",label=tags$b("Height   (m)"),value=values$Z)
  })
  output$uiD<-renderUI({
    numericInput("nD",label=tags$b("Diameter (m)"),value=values$D)
  })
  output$uiBn<-renderUI({
    sliderInput("nBn",label=tags$b("Blade Number (-)"),min=2,max=100,value=values$Bn,step=1)
  })
  output$uiH0<-renderUI({
    sliderInput("nH0",label=tags$b("Clearance Wall/Stirrer (m)"),min=0.001,max=0.1,value=values$H0,step=0.001)
  })
  output$uiRot<-renderUI({
    sliderInput("nRot",label=tags$b("Rotational Speed (rpm)"),min=100,max=1000,value=values$Rot,step=10)
  })
  output$uiP<-renderUI({
    numericInput("nP",label=tags$b("Pressure (bar)"),value=values$P,step=0.1)
  })
  output$uiTH<-renderUI({
    numericInput("nTH",label=tags$b("Wall Temperature (C)"),value=values$TH,step=1)
  })
  output$uiLFw<-renderUI({
    numericInput("nLFw",label=tags$b("Flow Rate   (Kg/h)"),value=values$LFw)
  })
  output$uiKc<-renderUI({
    numericInput("nhe",label = tags$b("Heat Transfer Efficiency (-)"),value=values$he)
  })
  output$uiRo<-renderUI({
    numericInput("nRo",label = tags$b("Density (kg/m3)"),value=values$Ro)
  })
  output$uiMu<-renderUI({
    numericInput("nMu",label = tags$b("Viscosity (Pa s)"),value=values$Mu)
  })
  output$uiW0<-renderUI({
    sliderInput("nW0",label=tags$b("Weight Dry Content (-)"),min=0,max=1,value=values$W0)
  })
  output$uiWF<-renderUI({
    sliderInput("nWF",label=tags$b("Weight Fraction Solvent (-)"),min=0,max=1,value=values$WF)
  })
  output$Nonsoleventcombo<-renderUI({
    selectInput("NonSolventSelector",label=tags$b("Remaining Solvent"),choices=list("None"=1,"water"=2),selected=values$NSS)
  })
  output$soleventcombo<-renderUI({
    choi<-as.list(1:nrow(values$df))
    names(choi)<-values$df[,1]
    selectInput("SolventSelector",label=tags$b("Extracting Solvent"),choices=choi,selected=values$SS)
  })
  output$examplecombo<-renderUI({
    values$Z<-values$ini[1,'nZ']
    values$D<-values$ini[1,'nD']
    values$Bn<-values$ini[1,'nBn']
    values$H0<-values$ini[1,'nH0']
    values$Rot<-values$ini[1,'nRot']
    values$P<-values$ini[1,'nP']
    values$TH<-values$ini[1,'nTH']
    values$LFw<-values$ini[1,'nLFw']
    values$he<-values$ini[1,'nhe']
    values$Ro<-values$ini[1,'nRo']
    values$Mu<-values$ini[1,'nMu']
    values$W0<-values$ini[1,'nW0']
    values$WF<-values$ini[1,'nWF']
    values$NSS<-as.numeric(values$ini[1,'NonSolventSelector'])
    values$SS<-as.numeric(values$ini[1,'SolventSelector'])
    choi<-as.list(1:nrow(values$ini))
    names(choi)<-values$ini[,1]
    selectInput("example",label=tags$b("Prebuilt Examples"),choices=choi,selected=as.numeric(values$ini[1,1]))
  })
  make_data<-reactive({
    df<-values$df
    n1<-values$SS
    pars<-c(
      MW1=df$MW[n1],                #kg/kmol
      D1=df$D[n1],                  #m2/s
      DH1=df$DH[n1],                #kcal/kmol
      Kc1=df$Kc[n1]/4180,           #kcal/ms
      roA1=df$roA[n1],
      roB1=df$roB[n1],
      ron1=df$ron[n1],
      roTc1=df$roTc[n1],
      CPA1=df$CpA[n1],
      CPB1=df$CpB[n1],
      CPC1=df$CpC[n1],
      CPD1=df$CpD[n1],
      P011=df$P01[n1],
      P012=df$P02[n1],
      P013=df$P03[n1],
      P014=df$P04[n1],
      P015=df$P05[n1],
      r1=df$r[n1],
      q1=df$q[n1],
      a12=df$a12[n1],
      a21=df$a21[n1]
    )
    n2<-values$NSS
    if(n2!=1){ # case of binary mixture with water
      pars<-c(pars,
        MW2=df$MW[1],
        D2=df$D[1],
        DH2=df$DH[1],
        Kc2=df$Kc[1]/4180,
        roA2=df$roA[1],
        roB2=df$roB[1],
        ron2=df$ron[1],
        roTc2=df$roTc[1],
        CPA2=df$CpA[1],
        CPB2=df$CpB[1],
        CPC2=df$CpC[1],
        CPD2=df$CpD[1],
        P021=df$P01[1],
        P022=df$P02[1],
        P023=df$P03[1],
        P024=df$P04[1],
        P025=df$P05[1],
        r2=df$r[1],
        q2=df$q[1]
      )
    }else{ # case of single organic solvent
      pars<-c(pars,
        MW2=1,
        DH2=0,
        Kc2=0,
        CPA2=0,
        CPB2=0,
        CPC2=0,
        CPD2=0,
        P021=0,
        P022=0,
        P023=0,
        P024=0,
        P025=0,
        r2=0,
        q2=0
      )
    }
    return(pars)
  })
  observeEvent(input$nZ,{
    req(input$nZ)
    values$Z<-as.numeric(input$nZ)
  })  
  observeEvent(input$nD,{
    req(input$nD)
    values$D<-as.numeric(input$nD)
  })  
  observeEvent(input$nH0,{
    req(input$nH0)
    values$H0<-as.numeric(input$nH0)
  })  
  observeEvent(input$nRot,{
    req(input$nRot)
    values$Rot<-as.numeric(input$nRot)
  })  
  observeEvent(input$nP,{
    req(input$nP)
    values$P<-as.numeric(input$nP)
  })  
  observeEvent(input$nTH,{
    req(input$nTH)
    values$TH<-as.numeric(input$nTH)
  })  
  observeEvent(input$nLFw,{
    req(input$nLFw)
    values$LFw<-as.numeric(input$nLFw)
  })  
  observeEvent(input$nRo,{
    req(input$nRo)
    values$Ro<-as.numeric(input$nRo)
  })  
  observeEvent(input$nhe,{
    req(input$nhe)
    values$he<-as.numeric(input$nhe)
  })  
  observeEvent(input$nMu,{
    req(input$nMu)
    values$Mu<-as.numeric(input$nMu)
  })  
  observeEvent(input$nW0,{
    req(input$nW0)
    values$W0<-as.numeric(input$nW0)
  })  
  observeEvent(input$nWF,{
    req(input$nWF)
    values$WF<-as.numeric(input$nWF)
  })  
  observeEvent(input$NonSolventSelector,{
    req(input$NonSolventSelector)
    values$NSS<-as.numeric(input$NonSolventSelector)
  })  
  observeEvent(input$SolventSelector,{
    req(input$SolventSelector)
    values$SS<-as.numeric(input$SolventSelector)
  })  
  observeEvent(input$example,{
    req(input$example)
    values$Z<-values$ini[as.numeric(input$example),'nZ']
    values$D<-values$ini[as.numeric(input$example),'nD']
    values$Bn<-values$ini[as.numeric(input$example),'nBn']
    values$H0<-values$ini[as.numeric(input$example),'nH0']
    values$Rot<-values$ini[as.numeric(input$example),'nRot']
    values$P<-values$ini[as.numeric(input$example),'nP']
    values$TH<-values$ini[as.numeric(input$example),'nTH']
    values$LFw<-values$ini[as.numeric(input$example),'nLFw']
    values$Kc<-values$ini[as.numeric(input$example),'nhe']
    values$Ro<-values$ini[as.numeric(input$example),'nRo']
    values$Mu<-values$ini[as.numeric(input$example),'nMu']
    values$W0<-values$ini[as.numeric(input$example),'nW0']
    values$WF<-values$ini[as.numeric(input$example),'nWF']
    values$NSS<-as.numeric(values$ini[as.numeric(input$example),'NonSolventSelector'])
    values$SS<-as.numeric(values$ini[as.numeric(input$example),'SolventSelector'])
  })  
  TFsimulation<-reactive({
    pars<-make_data()
    pars<-c(pars,
      SS=values$SS,
      NSS=values$NSS,
      Ro=values$Ro,
      he=values$he,
      Mu=values$Mu,
      W0=values$W0,
      LFw=values$LFw/3600,
      WF=values$WF,
      TH=values$TH+273,
      P=values$P,
      nrot=values$Rot/60,
      Bn=values$Bn,
      H0=values$H0,
      Z=values$Z,
      D=values$D,
      C=pi*values$D)
    if(values$NSS!=1){
          pars<-with(as.list(pars),
                     c(pars,
                       xF=WF/MW1/(WF/MW1+(1-W0-WF)/MW2),
                       LF=LFw*(WF/MW1+(1-W0-WF)/MW2)
                     )
          )
    }else{
          pars<-with(as.list(pars),
                     c(pars,
                       xF=1,
                       LF=LFw*(WF/MW1)
                     )
          )
    }
    dz<-log(1:60)
    z<-rev(values$Z*(dz[60]-dz)/dz[60])
    y<-with(c(as.list(pars)),{# initial conditions
      c(  L=LF,               #kmol/s
          x1=xF,              #-
          Qt=0,               #kcal/s
          DIS=0,              #kg/s
          DIS1=0)             #-
    })
    values$pars<-pars
    out<-lsodar(
      y=y,
      times=z,
      fun=Differential_system,
      rootfunc=Constrains,
      rtol=1e-4,
      atol=rep(1e-6,length(y)),
      parms=pars,
      events=list(func=eventfun,root=TRUE,terminalroot=c(1,2))
    )
    if (!is.null(attr(out,"indroot"))){
      if(attr(out,"indroot")==1){
        nid<-dim(out)[1]
        if(nid<length(z)){
          for(i in (nid+1):length(z)){
              out<-rbind(out,out[nid,])
          }
        }
        out<-as.data.frame(out)
        out$time<-z
        shinyalert(title = "All solvent is removed before the end of the column", type = "warning")
      }
    }
    #print(out)
    list(df=as.data.frame(out))
  })
  output$draw<- renderImage({list(src = "draw.jpg",width=500, height=500)}, deleteFile = FALSE)
  output$loadvalue<-renderText({
    validate(need(nrow(values$DS)!=0,""))
    paste(as.character(dim(values$DS)[1]*dim(values$DS)[2])," Values Loaded")
  })
  observeEvent(input$file_load,{
    req(input$file_load)
    inFile <- input$file_load
    if (is.null(inFile))
      return(NULL)
    req(input$rowindex)
    if (input$rowindex) {
      df <- read.csv(inFile$datapath,header = input$header,sep = input$sep,
                     row.names=1,quote = input$quote, dec=input$dec)
    }else{
      df <- read.csv(inFile$datapath,header = input$header,sep = input$sep,
                     quote = input$quote, dec=input$dec)
    }
    values$DS<-df
  })
  output$MFplot<-renderPlot({
    df<-TFsimulation()$df
    plt<-ggplot(df,aes(x=time,y=L))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    if(input$MFbox==2)plt<-plt+ scale_x_log10() 
    plt<-plt+labs(title="Liquid Molar Flow vs. Lenght", x="Lenght (m)", y ="Total Molar Flow (mol/s)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"),legend.position="none")
    print(plt)
  })
  output$HFplot<-renderPlot({
    df<-TFsimulation()$df
    plt<-ggplot(df,aes(x=time,y=Qt))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    if(input$HFbox==2)plt<-plt+ scale_x_log10() 
    plt<-plt+labs(title="Cumulate Heat Flow vs. Lenght", x="Lenght (m)", y ="Cumulate Heat Flow (Kcal/s)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"),legend.position="none")
    print(plt)
  })
  output$Tplot<-renderPlot({
    df<-TFsimulation()$df
    plt<-ggplot(df,aes(x=time,y=Teb-273.15))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    if(input$Tbox==2)plt<-plt+ scale_x_log10() 
    plt<-plt+labs(title="Temperature vs. Lenght", x="Lenght (m)", y ="Temperature (C)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"),legend.position="none")
    print(plt)
  })
  output$Cplot<-renderPlot({
    df<-TFsimulation()$df
    nr<-nrow(df)
    dfp<-data.frame(time=c(df$time,df$time,df$time),value=c(df$x1,df$x1eq,df$y1eq),colr=c(rep("a",nr),rep("b",nr),rep("c",nr)))
    dfp$colr<-factor(dfp$colr)
    plt<-ggplot(dfp,aes(x=time,y=value,colour=colr))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    if(input$Cbox==2)plt<-plt+ scale_x_log10() 
    plt<-plt+scale_color_manual(name="Molar Fractions",labels = c("Liquid Solvent","Liquid Solvent Equilibrium","Vapor Solvent Equilibrium"),values=c("a"="red","b"="blue","c"="green"))
    plt<-plt+ theme(legend.position="right")
    plt<-plt+labs(title="Solvent Molar Fraction in Liquid vs. Lenght", x="Lenght (m)", y ="Molar Fractions (mol/mol)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", colour="brown"))
    print(plt)
  })
  output$MTplot<-renderPlot({
    df<-TFsimulation()$df
    nr<-nrow(df)
    dfp<-data.frame(time=c(df$time,df$time,df$time),value=c(df$Q,df$Nt*1e5,df$N1*1e5),colr=c(rep("a",nr),rep("b",nr),rep("c",nr)))
    dfp$colr<-factor(dfp$colr)
    plt<-ggplot(dfp,aes(x=time,y=value,colour=colr))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    plt<-plt+scale_y_continuous(sec.axis = sec_axis(~./1e5, name = "Cross Mass Flow (mol/s)"))
    if(input$MTbox==2)plt<-plt+ scale_x_log10() 
    plt<-plt+ theme(legend.position="right")
    plt<-plt+scale_color_manual(name="Flow",labels = c("Heat","Mass Global","Mass Solvent"),values=c("a"="red","b"="blue","c"="green"))
    plt<-plt+labs(title="Mass Transfer vs. Lenght", x="Lenght (m)", y ="Cross Heat Flow (kcal/s)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"))
    print(plt)
  })
  tabsum<-reactive({
    out<-TFsimulation()$df
    nr<-nrow(out)
    with(as.list(values$pars),{
          Ex1<-(LF*xF-out$L[nr]*out$x1[nr])/LF/xF*100
          Pw<-out$Qt[nr]/C/Z*4.183*10
          xfin<-out$x1[nr]
          Lfin<-out$L[nr]*(xfin*MW1+(1-xfin)*MW2)+LFw*W0
          ER<-(LFw-Lfin)/LFw*100
          Gw<-(LFw-Lfin)*3600
          W0fin<-LFw*W0/Lfin
          W1Fin<-xfin*MW1/(xfin*MW1+(1-xfin)*MW2)
          Tfin<-out$Teb[nr]-273.15
          T0<-out$Teb[1]-273.15
          if(NSS==1){
            W10<-1
            W20<-0
            W2Fin<-0
            Ex2<-0
          }else{
            W10<-WF
            W20<-(1-WF)
            W2Fin<-(1-xfin)*MW2/(xfin*MW1+(1-xfin)*MW2)
            Ex2<-(LF*(1-xF)-out$L[nr]*(1-out$x1[nr]))/LF/(1-xF)*100
          } 
          topic<-c("Global Feed (kg/h)",
                   "Vapor Feed (kg/h)",
                   "Dry In (% w/w)",
                   "Solvent in Liquid (% w/w)",
                   "Non Solvent in Liquid (% w/w)",
                   "Boiling Temperature (C)",
                   "Pressure (bar)",
                   "Total Evaporation (% w)",
                   "Efficiency Solvent (% mol)",
                   "Efficiency Non Solvent (% mol)",
                   "Thermal Power(kW/m2)")
          valueIn<-c(toString(format(LFw*3600,digits=4)),
                     toString(format(0,digits=4)),
                     toString(format(W0*100,digits=4)),
                     toString(format(W10*100,digits=4)),
                     toString(format(W20*100,digits=4)),
                     toString(format(T0,digits=4)),
                     toString(format(P,digits=4)),
                     toString(format(0,digits=4)),
                     toString(format(0,digit=4)),
                     toString(format(0,digit=4)),
                     toString(format(0,digit=4)))
          valueOut<-c(toString(format(Lfin*3600,digits=4)),
                     toString(format(Gw,digits=4)),
                     toString(format(W0fin*100,digits=4)),
                     toString(format(W1Fin*100,digits=4)),
                     toString(format(W2Fin*100,digits=4)),
                     toString(format(Tfin,digits=4)),
                     toString(format(P,digits=4)),
                     toString(format(ER,digits=4)),
                     toString(format(Ex1,digit=4)),
                     toString(format(Ex2,digit=4)),
                     toString(format(Pw,digit=4)))
          return(data.frame(Topic=topic,Input=valueIn,Output=valueOut))
    })
  })
  output$viewsum<-({
    DT::renderDataTable({DT::datatable(tabsum(),caption=tags$b("Process Conditions and Outputs"), 
                                       class = 'cell-border stripe',
                                       options = list(paging = FALSE,searching = FALSE),
                                       filter='none',rownames=FALSE)})
  })
  tabflu<-reactive({
    out<-TFsimulation()$df
    Qv<-mean(out$Qv)
    topic<-c(  "Volumetric Flow (m3/h)",
               "Volumetric Fillet Flow (m3/h)",
               "Volumetric Film Flow (m3/h)",
               "Fillet Rate (m/s)",
               "Film Rate (m/s)",
               "Theoretical Thickeness (m)",
               "Film Thickeness (m)",
               "Fillet Diameter (m)",
               "Film Reynolds Number (-)",
               "Global Heat Exchange Coefficient (W/m2K)",
               "Global Mass Exchange Coefficient (kg/hm2)")
      value<-c( toString(format(mean(out$Qv)*3600,digits=4)),
                toString(format(mean(out$Qf)*3600,digits=4)),
                toString(format(mean(out$Qm)*3600,digits=4)),
                toString(format(mean(out$wf),digits=4)),
                toString(format(mean(out$wm),digits=4)),
                toString(format(mean(out$Hv),digits=4)),
                toString(format(mean(out$Hm),digits=4)),
                toString(format(mean(out$Hf),digits=4)),
                toString(format(mean(out$ReL),digit=4)),
                toString(format(mean(out$h)*4180,digit=4)),
                toString(format(mean(out$kx)*3600,digit=4)))
      return(data.frame(Topic=topic,Value=value))
  })
  output$viewflu<-({
    DT::renderDataTable({DT::datatable(tabflu(),caption=tags$b("Flow,Thickness and Coefficient"), 
                                       class = 'cell-border stripe',
                                       options = list(paging = FALSE,searching = FALSE),
                                       filter='none',rownames=FALSE)})
  })
  observeEvent(input$file_load,{
    req(input$file_load,input$rowindex)
    inFile <- input$file_load
    if (input$rowindex) {
      values$dfopt<- read.csv(inFile$datapath,header = input$header,sep = input$sep,
                     row.names=1,quote = input$quote, dec=input$dec)
    }else{
      values$dfopt<- read.csv(inFile$datapath,header = input$header,sep = input$sep,
                     quote = input$quote, dec=input$dec)
    }
  })
  output$xfplot<-renderPlot({
    RMSE<-rep(0,input$eslider[2]-input$eslider[1]+1)
    np<-nrow(values$dfopt)
    vHE<-input$eslider[1]:input$eslider[2]
    ne<-length(vHE)
    j<-0
    HEmin<-vHE[1]
    RMSEmin<-100
    withProgress(message = 'Making plot',value=0,min=0,max=ne,{
            for(j in 1:ne){
              REc<-rep(0,np)
              for (nt in 1:np){
                values$SS<-as.numeric(values$dfopt[nt,'SolventSelector'])
                values$NSS<-as.numeric(values$dfopt[nt,'NonSolventSelector'])
                pars<-make_data()
                pars<-c(pars,
                        SS=values$SS,
                        NSS=values$NSS,
                        he=vHE[j]/100,
                        Ro=values$dfopt[nt,'nRo'],
                        Mu=values$dfopt[nt,'nMu'],
                        W0=values$dfopt[nt,'nW0'],
                        LFw=values$dfopt[nt,'nLFw']/3600,
                        WF=values$dfopt[nt,'nWF'],
                        TH=values$dfopt[nt,'nTH']+273,
                        P=values$dfopt[nt,'nP'],
                        nrot=values$dfopt[nt,'nRot']/60,
                        Bn=values$dfopt[nt,'nBn'],
                        H0=values$dfopt[nt,'nH0'],
                        Z=values$dfopt[nt,'nZ'],
                        D=values$dfopt[nt,'nD'],
                        C=pi*values$dfopt[nt,'nD'])
                if(values$NSS!=1){
                  pars<-with(as.list(pars),
                             c(pars,
                               xF=WF/MW1/(WF/MW1+(1-W0-WF)/MW2),
                               LF=LFw*(WF/MW1+(1-W0-WF)/MW2)
                             )
                  )
                }else{
                  pars<-with(as.list(pars),
                             c(pars,
                               xF=1,
                               LF=LFw*(WF/MW1)
                             )
                  )
                }
                dz<-log(1:60)
                z<-rev(values$dfopt[nt,'nZ']*(dz[60]-dz)/dz[60])
                y<-with(c(as.list(pars)),{# initial conditions
                  c(  L=LF,               #kmol/s
                      x1=xF,              #-
                      Qt=0,               #kcal/s
                      DIS=0,              #kg/s
                      DIS1=0)             #-
                })
                values$pars<-pars
                out<-lsodar(
                  y=y,
                  times=z,
                  fun=Differential_system,
                  rootfunc=Constrains,
                  rtol=1e-4,
                  atol=rep(1e-6,length(y)),
                  parms=pars,
                  events=list(func=eventfun,root=TRUE,terminalroot=c(1,2))
                )
                if (!is.null(attr(out,"indroot"))){
                  if(attr(out,"indroot")==1){
                    nid<-dim(out)[1]
                    if(nid<length(z)){
                      for(i in (nid+1):length(z)){
                        out<-rbind(out,out[nid,])
                      }
                    }
                    out<-as.data.frame(out)
                    out$time<-z
                  }
                }
                out<-as.data.frame(out)
                REc[nt]<-(out$L[nrow(out)]/out$L[1]-values$dfopt$nWE[nt])^2
              }
              RMSE[j]<-sum(REc)
              if(RMSE[j]<RMSEmin){
                RMSEmin<-RMSE[j]
                HEmin<-vHE[j]
              }
              incProgress(1, detail = paste("Doing part", j))
            }
    })
    values$HEopt<-HEmin
    dfxf<-data.frame(HE=input$eslider[1]:input$eslider[2],RMSE=RMSE)
    plt<-ggplot(dfxf,aes(x=HE,y=RMSE))+theme_bw()
    plt<-plt+geom_line()+geom_point()
    plt<-plt+labs(title=paste("RMSE vs. Parameter : minimum in :",format(HEmin,digits=4),sep=' '), x="Efficiency (-)", y ="RMSE (-)")
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"),legend.position="none")
    print(plt)
  })
  respfunc<-function(he){
    np<-nrow(values$dfopt)
    WEc<-rep(0,np)
    for (nt in 1:np){
      values$SS<-as.numeric(values$dfopt[nt,'SolventSelector'])
      values$NSS<-as.numeric(values$dfopt[nt,'NonSolventSelector'])
      pars<-make_data()
      pars<-c(pars,
              SS=values$SS,
              NSS=values$NSS,
              he1=he[1],
              he2=he[2],
              he3=he[3],
              he4=he[4],
              he5=he[5],
              he6=he[6],
              Ro=values$dfopt[nt,'nRo'],
              Mu=values$dfopt[nt,'nMu'],
              W0=values$dfopt[nt,'nW0'],
              LFw=values$dfopt[nt,'nLFw']/3600,
              WF=values$dfopt[nt,'nWF'],
              TH=values$dfopt[nt,'nTH']+273,
              P=values$dfopt[nt,'nP'],
              nrot=values$dfopt[nt,'nRot']/60,
              Bn=values$dfopt[nt,'nBn'],
              H0=values$dfopt[nt,'nH0'],
              Z=values$dfopt[nt,'nZ'],
              D=values$dfopt[nt,'nD'],
              C=pi*values$dfopt[nt,'nD'])
      if(values$NSS!=1){
        pars<-with(as.list(pars),
                   c(pars,
                     xF=WF/MW1/(WF/MW1+(1-W0-WF)/MW2),
                     LF=LFw*(WF/MW1+(1-W0-WF)/MW2)
                   )
        )
      }else{
        pars<-with(as.list(pars),
                   c(pars,
                     xF=1,
                     LF=LFw*(WF/MW1)
                   )
        )
      }
      dz<-log(1:60)
      z<-rev(values$dfopt[nt,'nZ']*(dz[60]-dz)/dz[60])
      y<-with(c(as.list(pars)),{# initial conditions
        c(  L=LF,               #kmol/s
            x1=xF,              #-
            Qt=0,               #kcal/s
            DIS=0,              #kg/s
            DIS1=0)             #-
      })
      values$pars<-pars
      out<-lsodar(
        y=y,
        times=z,
        fun=Differential_system,
        rootfunc=Constrains,
        rtol=1e-4,
        atol=rep(1e-6,length(y)),
        parms=pars,
        events=list(func=eventfun,root=TRUE,terminalroot=c(1,2))
      )
      if (!is.null(attr(out,"indroot"))){
        if(attr(out,"indroot")==1){
          nid<-dim(out)[1]
          if(nid<length(z)){
            for(i in (nid+1):length(z)){
              out<-rbind(out,out[nid,])
            }
          }
          out<-as.data.frame(out)
          out$time<-z
        }
      }
      out<-as.data.frame(out)
      WEc[nt]<-out$L[nrow(out)]/out$L[1]
    }
    #print(WEc)
  return(WEc) 
  }
  optfunc<-function(he){
    WEm<-values$dfopt$nWE
    WEc<-respfunc(he)
    RMSE<-sum((WEc-WEm)^2)
    print(c(he,RMSE))
    return(RMSE)
  }
  output$xfres<-renderPlot({
    #h<-he1*nrot^he2*TH^he3*LFw^he4    
    # he1<-1e-4
    # he2<-0.5
    # he3<-1
    # he4<-0.5
    # he1<-0.0000117331
    # he2<-0.0571532473
    # he3<-2.2629273020
    # he4<-0.6847104698
    # he1<-5.583121e-04
    # he2<-1.363871e-03
    # he3<-1.510147e-04
    # he4<--6.320389e-04
    # he5<-8.433021e-06
    # he6<--9.520139e-04
    he1<-5.583121e-04
    he2<-1.363871e-03
    he3<-1.510147e-04
    opt_out<-optim_nm(fun=optfunc,k=3,start=c(he1,he2,he3),trace=TRUE)
    print(opt_out)
    df<-as.data.frame(opt_out$trace)
    WEc<-respfunc(opt_out$par)
    WEm<-values$dfopt$nWE
    df<-data.frame(WEc=WEc,WEm=WEm)
    plt<-ggplot(df,aes(x=WEm,y=WEc))+theme_bw()
    plt<-plt+geom_point()+xlim(min(df),max(df))+ylim(min(df),max(df))
    plt<-plt+labs(title="Calculated vs. Measured",x="X measured", y ="X calulated")
    plt<-plt+ geom_abline(intercept=0,slope=1)
    plt<-plt+theme(plot.title=element_text(face="bold",size="14", color="brown"),legend.position="none")
    print(plt)
  })
  output$hlp_geometry<-renderUI({includeHTML("hlp_geometry.html")})
  output$hlp_process<-renderUI({includeHTML("hlp_process.html")})
  output$hlp_expdata<-renderUI({includeHTML("hlp_expdata.html")})
  output$hlp_feedglob<-renderUI({includeHTML("hlp_feedglob.html")})
  output$hlp_feedcomp<-renderUI({includeHTML("hlp_feedcomp.html")})
  output$hlp_MFsimul<-renderUI({includeHTML("hlp_MFsimul.html")})
  output$hlp_Fsimul<-renderUI({includeHTML("hlp_Fsimul.html")})
  output$hlp_summary<-renderUI({includeHTML("hlp_summary.html")})
  output$hlp_xfopt<-renderUI({includeHTML("hlp_xfopt.html")})
  output$hlp_xfres<-renderUI({includeHTML("hlp_xfres.html")})
} # end of server
