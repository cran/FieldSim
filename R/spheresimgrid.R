### spheresimgrid.R  (2006-10)
###
###    Random field simulation on the sphere by the FieldSim method
###
### Copyright 2006-2070 Alexandre Brouste and Sophie Lambert-Lacroix
###
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

spheresimgrid <- function(R,Ne=100,Nr=100,nbNeighbor=4,S){

##    INPUT VARIABLES
#########################
##  R   : R function
##      covariance function
##  Ne  : positive integer
##      Number of simulation points associated with the uniform distributed discretization 
##      of the sphere for the first step of the algorithm (Accurate simulation step)
##  Nr  : positive integer
##      Number of simulation points associated with the uniform distributed discretization 
##      of the sphere for the second step of the algorithm (Refined simulation step).
##  nbNeighbor  : positive integer
##      number of neighbors to use in the second step 
##      of the algorithm. It must be between 1 and 32.
##  S  : matrix of size 3xN
##      S[,n] gives the sphere point (S[1,n],S[2,n],S[3,n]) where the field must be 
##      simulated


##    OUTPUT VARIABLES
##########################
##  S  : matrix of size 3xN
##      S[,n] gives the sphere point (S[1,n],S[2,n],S[3,n]) where the field has been 
##      simulated
##
##  Wg  : vector of length N; Wg[n] gives the value of the process
##  at point (S[1,n],S[2,n],S[3,n]).
##      
##  time: CPU time

##  TEST ON INPUT VARIABLES
##############################

if (is.function(R)==FALSE){
 stop("Message from spheresimgrid.R: R is not of valid type")}
 
if ((is.numeric(Ne)==FALSE)||(round(Ne)-Ne!=0)||(Ne<0)){
 stop("Message from spheresimgrid.R: Ne is not of valid type")}
 
if ((is.numeric(Nr)==FALSE)||(round(Nr)-Nr!=0)||(Nr<0)){
 stop("Message from spheresimgrid.R: Nr is not of valid type")}
 
if ((is.numeric(nbNeighbor)==FALSE)||(round(nbNeighbor)-nbNeighbor!=0)||(nbNeighbor<1)||(nbNeighbor>32)){
 stop("Message from spheresimgrid.R: nbNeighbor is not of valid type")}

if (nbNeighbor>Ne+1){
 stop("Message from spheresimgrid.R: nbNeighbor is too large or Ne is to small")}

if((is.matrix(S)==FALSE)||(nrow(S)!=3)||(apply(S^2,2,sum)>rep(1+10^(-14),ncol(S)))||(apply(S^2,2,sum)<rep(1-10^(-14),ncol(S)))){
stop("Message from spheresimgrid.R: S is not of valid type")}

#-----------------------------------------------------
#------Accurate simulation step for j<Ne    ------
#-----------------------------------------------------

time <- proc.time()

#Simulation of Ne+Nr points of discretization on the sphere
Theta <- runif(Ne+Nr, min=0, max=2*pi)
Phi <- acos(1-2*runif(Ne+Nr, min=0, max=1))
X <- cos(Theta)*sin(Phi)
Y <- sin(Theta)*sin(Phi)
Z <- cos(Phi)
W <- rep(0,length(Z))




J<-Ne #discretization level
#******initialization*******

    ts<-1                     #simulation time
    tr<-c(X[ts],Y[ts],Z[ts])  #spatial abcisse
    sigma2f<-R(cbind(tr,tr)) 
    sigma2bold<-sigma2f 
    W[ts]<-sqrt(sigma2f)*rnorm(1)
  
#******next steps***********
for (ts in 2:J){

    tr_x<-X[ts]
    tr_y<-Y[ts]
    tr_z<-Z[ts]
    tr<-c(tr_x,tr_y,tr_z)
   
    #Compute W(ts)

    #///////////////////////////////////////////
    sigma2f<-R(cbind(tr,tr))        #variance //
    Af<-rep(0,(ts-1))   #filtre de taille ts-1//
    #///////////////////////////////////////////

    #******init****************
    sigma2bnew<-sigma2f
    tabR<-sigma2f
    Abnew<-matrix(0,(ts-1),(ts-1))
    index<-NULL
    indextemp<-index

    for (k in 1:(ts-1)){
       indextemp<-c(k,indextemp)
       #*****Compute the filters Af_k ***
       t<-R(cbind(tr,c(X[(ts-k)],Y[(ts-k)],Z[(ts-k)])))  
       tabR<-c(t,tabR) #taille k+1
       s<-sigma2bold[k]
       #//// Compute beta(ts,ts-k)

       if (k==1)
       stdaux <- sqrt(sigma2f*s)
       if (stdaux==0)
        {beta <- 0}
       if (stdaux!=0)
        {beta <- t/stdaux}
        
       if (k!=1){
       b <- Abold[k-1,1:(k-1)]
       stdaux <- sqrt(sigma2f*s)
       if (stdaux==0)
        {beta <- 0}
       if (stdaux!=0)
        {
         beta <- (t+ sum(b*tabR[2:(length(tabR)-1)]))/stdaux}
       }
       
       #
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

     #********end of computation of filter Af ***
     }

     #stockage
     sigma2bold<-sigma2bnew
     Abold<-Abnew


  #/////////Simul W(ts) ///////////////////
    
  
    W[ts]<-sqrt(sigma2f)*rnorm(1) - sum(Af*W[indextemp])            
  }
#----------------------------------------------------
#------end of the accurate simulation step ----------
#----------------------------------------------------


#----------------------------------------------------
#----       approximation step     ------------------
#----------------------------------------------------




for (ts in (J+1):(Nr+Ne)){


    #////////compute the indexes of neighbors //////////
    u<-matrix(c(X[ts],Y[ts],Z[ts]),nrow=1,ncol=3)%*%rbind(X[1:(ts-1)],Y[1:(ts-1)],Z[1:(ts-1)])
    u[u>1] <- 1
    u[u<(-1)] <- -1
    D <- acos(u)
    Voisins <- sort(D,index=TRUE)$ix[1:nbNeighbor]
    
    
    #////////compute filter and variance //////
    xx<-rbind(c(X[Voisins],X[ts]),c(Y[Voisins],Y[ts]),c(Z[Voisins],Z[ts]))
    out<-.Call("vf3",as.double(xx),body(R),new.env())
    var<-out[1] 
    filtre<-out[2:length(out)]

    #////////simulation of W(ts) /////////////////////
    W[ts]<-sqrt(var)*rnorm(1) - sum(filtre[1:length(Voisins)]*W[Voisins])    
}

#----------------------------------------------------
#------end of the approximation simulation step -----
#----------------------------------------------------

#----------------------------------------------------
#----       final approximation step     -----------
#----------------------------------------------------
N <- ncol(S)
X <- c(X,S[1,])
Y <- c(Y,S[2,])
Z <- c(Z,S[3,])
Wg <- numeric(0)
for (ts in (Nr+Ne+1):(Nr+Ne+N)){

    #////////compute the indexes of neighbors //////////
    u<-matrix(c(X[ts],Y[ts],Z[ts]),nrow=1,ncol=3)%*%rbind(X[1:(ts-1)],Y[1:(ts-1)],Z[1:(ts-1)])
    u[u>1] <- 1
    u[u<(-1)] <- -1
    D <- acos(u)
    Voisins <- sort(D,index=TRUE)$ix[1:nbNeighbor]
    
    
    #////////compute filter and variance //////
    xx<-rbind(c(X[Voisins],X[ts]),c(Y[Voisins],Y[ts]),c(Z[Voisins],Z[ts]))
    out<-.Call("vf3",as.double(xx),body(R),new.env())
    var<-out[1] 
    filtre<-out[2:length(out)]

    #////////simulation of W(ts) /////////////////////
    W[ts]<-sqrt(var)*rnorm(1) - sum(filtre[1:length(Voisins)]*W[Voisins])    
    Wg <- c(Wg,W[ts])
}

#--------------------------------------------------------------
#------end of the final approximation simulation step --------
#--------------------------------------------------------------

time <- proc.time() - time

return(list(S=S,Wg=Wg,time=time[[3]]))
}
