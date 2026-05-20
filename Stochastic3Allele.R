## relatively simple, stochastic 3 allele model with vs without NFDS
## f = initial frequency of alleles
## w = relative fitness (if nfds != this if relative fitness at p = 0.5)
## nfds = effect of frequency on relative fitness
## N = number of gens to simulate
## odm = fitness benefit for melanic hets
## Ne = 110
## probNFDS = 1, probability that NFDS occurs
sim3a<-function(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=1,wm=1,nfdss=0,nfdsg=0,nfdsm=0,N=500,odm=0,Ne=110,probNFDS=1){
    P<-matrix(NA,nrow=N,ncol=3) ## s, g, m freqs
    P[1,]<-c(fs,fg,fm)
    
    
    M<-P ## morph frequencies, dominance g > s > m
    W<-P  ## relative fitness values, s, g, m

    for(j in 2:N){
    	ppp<-P[j-1,]
    	yyy<-rmultinom(n=1,size=Ne*2,prob=ppp)
        fs<-yyy[1]/(2*Ne)
        fg<-yyy[2]/(2*Ne)
        fm<-yyy[3]/(2*Ne)
        M[j-1,1]<-fs^2 + 2*fs*fm
        M[j-1,2]<-fg^2 + 2*fg*fm + 2*fs*fg
        M[j-1,3]<-fm^2

        w<-rep(NA,3)
        snfds<-rbinom(1,1,prob=probNFDS)
        w[1]<-ws + nfdss *(M[j-1,1]-.5) * snfds   
        w[2]<-wg + nfdsg *(M[j-1,2]-.5) * snfds    
        w[3]<-wm + nfdsm *(M[j-1,3]-.5) * snfds           
        w[w<0]<-0
        w<-w/(w[1]*M[j-1,1] + w[2]*M[j-1,2] + w[3]*M[j-1,3])
        W[j-1,]<-w                
        fsp<-fs^2*w[1] + fs*fg*w[2] + fs*fm*(w[1] + odm)
        fgp<-fs*fg*w[2] + fg^2*w[2] + fg*fm*(w[2] + odm)
        fmp<-fs*fm*(w[1] + odm) + fg*fm*(w[2] + odm) + fm^2*w[3]
        sump<-fsp+fgp+fmp
        P[j,1]<-fsp/sump
        P[j,2]<-fgp/sump
        P[j,3]<-fmp/sump
     }
     out<-cbind(P,M,W)
     return(out)
     }
 

plotsim<-function(o=NULL,x=0,tit=NULL,w=c(1,1,1),cl=1.2,ct=1.1){##x to 3 for morphs
    plot(o[,1+x],ylim=c(0,1),type='l',col="cadetblue",xlab="Generation",ylab="Frequency",cex.lab=cl)
    lines(o[,2+x],col="forestgreen")
    lines(o[,3+x],col="brown")
    legend(-.03,1.03,w,fill=c("cadetblue","forestgreen","brown"),bty='n',title="fitness",ncol=3,cex=ct)
    title(main=tit)
}

pdf("SSim3A_freqs.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,4.5,2.5,1))
 
oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=0,nfdsg=0,nfdsm=0)
plotsim(oo,x=0,tit="No NFDS",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.5,nfdsg=-1.5)
plotsim(oo,x=0,tit="NFDS G+S (-1.5)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,N=25000)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=2/3,fg=1/3,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=1/2,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/4,fg=3/4,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=0,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=0,tit="NFDS G+S+M (-1.1)")

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=0,tit="NFDS G+S+M (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.7,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=0,tit="NFDS G+S+M (-1.1)",w=c(1,.7,.95))

dev.off()

pdf("SSim3A_morphs.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,4.5,2.5,1))
 
oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=0,nfdsg=0,nfdsm=0)
plotsim(oo,x=3,tit="No NFDS",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.5,nfdsg=-1.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.5)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,N=25000)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=2/3,fg=1/3,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=1/2,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/4,fg=3/4,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1)
plotsim(oo,x=3,tit="NFDS G+S (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1)")

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.7,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1)",w=c(1,.7,.95))

dev.off()

pdf("SSim3A_morphs_od.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,4.5,2.5,1))
 
oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=0,nfdsg=0,nfdsm=0,odm=0.05)
plotsim(oo,x=3,tit="OD only (.05)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.5,nfdsg=-1.5,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.5) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,N=25000,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=2/3,fg=1/3,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=1/2,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/4,fg=3/4,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.7,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.7,.95))

dev.off()

pdf("SSim3A_morphs_od_ceanothus.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,4.5,2.5,1))
 
oo<-sim3a(ws=0.9,wg=1,wm=.95,nfdss=0,nfdsg=0,nfdsm=0,odm=0.05)
plotsim(oo,x=3,tit="OD only (.05)",w=c(0.9,1,.95))

oo<-sim3a(ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(ws=0.9,wg=1,wm=.95,nfdss=-1.5,nfdsg=-1.5,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.5) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,N=25000,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=2/3,fg=1/3,fm=0,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/2,fg=1/2,fm=0,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/4,fg=3/4,fm=0,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=0.9,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(0.9,1,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=0.7,wg=1,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(.7,1,.95))

dev.off()

pdf("SSim3A_morphs_od_nfds50.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,4.5,2.5,1))
 
oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=0,nfdsg=0,nfdsm=0,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="OD only (.05)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(ws=1,wg=0.9,wm=.95,nfdss=-1.5,nfdsg=-1.5,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.5) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=4/6,fg=1/6,fm=1/6,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,N=25000,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=2/3,fg=1/3,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=1/2,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/4,fg=3/4,fm=0,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/2,fg=0,fm=1/2,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.9,.95))

oo<-sim3a(fs=1/3,fg=1/3,fm=1/3,ws=1,wg=0.7,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=0.05,probNFDS=.5)
plotsim(oo,x=3,tit="NFDS G+S+M (-1.1) + OD (.05)",w=c(1,.7,.95))

dev.off()
