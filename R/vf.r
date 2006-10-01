VF <-function(x,R){
z <- dim(x)[1]
ts <- dim(x)[2]
u<-x[1:2,1]
sigma2f<-R(cbind(u,u))

sigma2bold<-sigma2f

    for (kk in 2:ts){ # pour tous les deltas du tableau

         u<-x[1:2,kk]
         #//////////////////////////////////////////////
         sigma2f<-R(cbind(u,u)) 
  #variance                  //
         Af<-rep(0,(kk-1))     #filtre de taille k-1//
         #//////////////////////////////////////////////

         #******init****************

         sigma2bnew<-sigma2f
         tabR<-sigma2f
         Abnew<-matrix(0,(kk-1),(kk-1))
         index<-NULL

         for (k in 1:(kk-1)){

           #*****calcul des filtres Af_k ***
           u<-x[1:2,kk]
           v<-x[1:2,(kk-k)]
           t<-R(cbind(u,v))
           tabR<-c(t,tabR) #taille k+1
           s<-sigma2bold[k]
#print(s)
   

           #//// calcul de betadelta(u,v)

           if (k==1){
           stdaux <- sqrt(sigma2f*s)
           if (stdaux ==0)
            {beta<-0}
           if (stdaux !=0)
            {beta<-t/stdaux}
           }
           if (k!=1){
           b<-Abold[(k-1),1:(k-1)]
           stdaux <- sqrt(sigma2f*s)
           if (stdaux ==0)
            {beta<-0}
           if (stdaux !=0)
            {beta<-(t+ sum(b*tabR[2:(length(tabR)-1)]))/stdaux}
           }

           if (abs(beta)>1){
            stop("message from FieldSim: R is not positive defined")
           }

           if (stdaux ==0){
            Af[k]<-0
            Abnew[k,k]<-0}

           if (stdaux !=0){
            r<-sqrt(sigma2f/s)
            Af[k]<--beta*r
            Abnew[k,k]<--beta/r}

           if (k>1){
           index<-c(k-1,index)
           a<-Af[1:(k-1)]
           Af[1:(k-1)]<-a+Af[k]*b[index]
           Abnew[k,1:(k-1)]<-b+Abnew[k,k]*a[index]
           }

           #update variance
           sigma2f<-(1-beta^2)*sigma2f
           sigma2bnew<-c(sigma2bnew,(1-beta^2)*s)

           #********fin calcul filtre Af ***
         } #end for k

         #stockage
         sigma2bold<-sigma2bnew
         Abold<-Abnew
         } #end for kk


return(list(variance=sigma2f,filtre=rev(Af[1:(ts-1)])))}
