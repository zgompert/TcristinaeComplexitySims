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
 
so1<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=1,wg=.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,odm=.05,N=100000,Ne=110)
    so1[i]<-sum(o[100000,1:3] > 0.0001)
    }
    
so2<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=1,wg=.9,wm=.95,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=.05,N=100000,Ne=110)
    so2[i]<-sum(o[100000,1:3] > 0.0001)
    }
    
so3<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=1,wg=.7,wm=.8,nfdss=-1.1,nfdsg=-1.1,nfdsm=-1.1,odm=.05,N=100000,Ne=110)
    so3[i]<-sum(o[100000,1:3] > 0.0001)
    }        

so4<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=.4,wg=.2,wm=.2,nfdss=-.5,nfdsg=-.1,nfdsm=0,odm=.025,N=100000,Ne=110)
    so4[i]<-sum(o[100000,1:3] > 0.0001)
    }
    
so5<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=.4,wg=.2,wm=.2,nfdss=-.5,nfdsg=-.1,nfdsm=0,odm=.025,N=100000,Ne=1000)
    so5[i]<-sum(o[100000,1:3] > 0.0001)
    }    

so6<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=.4,wg=.2,wm=.2,nfdss=-.5,nfdsg=-.1,nfdsm=-.1,odm=.025,N=100000,Ne=110)
    so6[i]<-sum(o[100000,1:3] > 0.0001)
    }    
         
         
so7<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=.4,wg=.2,wm=.2,nfdss=-.5,nfdsg=-.1,nfdsm=-.1,odm=.025,N=100000,Ne=10000)
    so7[i]<-sum(o[100000,1:3] > 0.0001)
    }    
        
so7<-rep(NA,20)
for(i in 1:20){
    o<-sim3a(ws=.4,wg=.2,wm=.15,nfdss=-.5,nfdsg=-.1,nfdsm=-.1,odm=.05,N=100000,Ne=110)
    so7[i]<-sum(o[100000,1:3] > 0.0001)
    }    
                  

          
