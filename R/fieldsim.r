### fieldsim.R  (2006-10)
###
###    Random field simulation by the method FieldSim method
###
### Copyright 2006-10 Alexandre Brouste and Sophie Lambert-Lacroix
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

fieldsim <- function(R,Elevel=1,Rlevel=5,nbNeighbor=4){

##    INPUT VARIABLES
#########################
##  R   : R function
##      covariance function
##  Elevel  : positive integer
##      Level associated with the regular space discretization 
##      for the first step of the algorithm (Accurate simulation step)
##      Form of regular space discretization : [[0:2^(Elevel)]/2^(Elevel)]^2
##  Rlevel  : positive integer
##      Level associated with the regular space discretization 
##      for the second step of the algorithm (Refined simulation step).
##      Form of regular space discretization : [[0:2^(Elevel+Rlevel)]/2^(Elevel+Rlevel)]^2
##  nbNeighbor  : positive integer
##      number of neighbors to use in the second step 
##      of the algorithm. It must be between 1 and 32.


##    OUTPUT VARIABLES
##########################
## Zrow: vector of length 2^(Elevel+Rlevel)+1
##      Zrow = 0:2^(Elevel+Rlevel)]/2^(Elevel+Rlevel)
## Zcol: vector of length 2^(Elevel+Rlevel)+1
##      Zcol = 0:2^(Elevel+Rlevel)]/2^(Elevel+Rlevel)
## Z: matrix of size (2^(Elevel+Rlevel)+1)x(2^(Elevel+Rlevel)+1)
##      Z[i,j] gives the value of the process at point (Zrow[i],Zcol[j])
## time: CPU time

##  TEST ON INPUT VARIABLES
##############################

if (is.function(R)==FALSE){
 stop("Message from fieldsim.R: R is not of valid type")}
 
if ((is.numeric(Elevel)==FALSE)||(round(Elevel)-Elevel!=0)||(Elevel<0)){
 stop("Message from fieldsim.R: Elevel is not of valid type")}
 
if ((is.numeric(Rlevel)==FALSE)||(round(Rlevel)-Rlevel!=0)||(Rlevel<0)){
 stop("Message from fieldsim.R: Rlevel is not of valid type")}
 
if ((is.numeric(nbNeighbor)==FALSE)||(round(nbNeighbor)-nbNeighbor!=0)||(nbNeighbor<1)||(nbNeighbor>32)){
 stop("Message from fieldsim.R: nbNeighbor is not of valid type")}

if (nbNeighbor>(2^(Elevel)+1)^2){
 stop("Message from fieldsim.R: nbNeighbor is too large or Elevel is to small")}



#-----------------------------------------------------
#------Accurate simulation step for j<Elevel    ------
#-----------------------------------------------------

time <- proc.time()
J<-Elevel #discretization level
R0 <- R(cbind(0,0,0,0))

#******initialization*******
for (l in 0:2^J){
 for (m in 0:2^J){
 
  if ((R0!=0)&(m==0)&(l==0)){
    ts<-1                     #simulation time
    tr<-c(0,0)                #spatial abcisse
    sigma2f<-R(cbind(tr,tr)) 
    sigma2bold<-sigma2f 
    x<-sqrt(sigma2f)*rnorm(1)
    X<-c(tr,x)
  }
 
  if ((l!=0) | (m!=0)){
   if((l==0) & (m==1) & (R0==0)){  
    ts<-1                     #simulation time
    tr<-c(0,1/2^J)            #spatial abcisse
    sigma2f<-R(cbind(tr,tr)) 
    sigma2bold<-sigma2f 
    x<-sqrt(sigma2f)*rnorm(1)
    X<-c(tr,x)
   }
   
   if (((R0==0)&((l!=0)|(m!=1)))|(R0!=0)){
    #for ts = 2..2^J

    ts<-ts+1
    tr_x<-l/2^(J)
    tr_y<-m/2^(J)
    tr<-c(tr_x,tr_y)
   
    #Compute x(ts)

    #///////////////////////////////////////////
    sigma2f<-R(cbind(tr,tr))        #variance     //
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
       if (is.vector(X)==TRUE)
         t<-R(cbind(tr,X[1:2]))
       if (is.matrix(X)==TRUE)
         t<-R(cbind(tr,X[1:2,(ts-k)]))  
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


  #/////////Simul X(tn) ///////////////////
    
  if (is.vector(X)==TRUE)
    x<-sqrt(sigma2f)*rnorm(1) - sum(Af*X[3])
  if (is.matrix(X)==TRUE)  
    x<-sqrt(sigma2f)*rnorm(1) - sum(Af*X[3,indextemp])            
  X<-cbind(X,c(tr,x))
  }
 }
}
}

if (R0==0){
 X<-cbind(c(0,0,0),X)} 
ts<-ts+1 

Z<-matrix(X[3,],2^J+1,2^J+1,byrow=TRUE)

#----------------------------------------------------
#------end of the accurate simulation step ----------
#----------------------------------------------------

Zold<-Z

#----------------------------------------------------
#----       approximation step     ------------------
#----------------------------------------------------


niveau <- J+1

if (nbNeighbor<=4)
    d <- 1
if ((nbNeighbor>4)&(nbNeighbor<=8))
    d <- 3
if ((nbNeighbor>8)&(nbNeighbor<=16))
    d <- 5
if ((nbNeighbor>16)&(nbNeighbor<=32))
    d <- 7


while (niveau<=(Elevel+Rlevel)){
Y <- matrix(0,2^(niveau)+1,2^(niveau)+1)

#------Build the matrix Y-------------------------------

for (l in 1:(2^(niveau)+1)){ #columns
for (m in 1:(2^(niveau)+1)){ #rows
if (((m/2-floor(m/2))!=0) & ((l/2-floor(l/2))!=0))
Y[m,l]<-Z[((m-1)/2+1),((l-1)/2+1)] 
}}

for (m in 1:2^(niveau-1)){ 
for (l in 1:2^(niveau-1)) {
   
    ###################################
    #Simulation for the midle point 
    ###################################
    pc_x<-2*l
    pc_y<-2*m
    tr_x<-(pc_x-1)/2^(niveau)
    tr_y<-(pc_y-1)/2^(niveau)
    
    tc<-c(tr_x,tr_y)

    
    #Compute Yaux the submatrix of Y containing the midle point neighbors
    Indexrow <- max(1,pc_x-d):min(2^(niveau)+1,pc_x+d)
    Indexcol <- max(1,pc_y-d):min(2^(niveau)+1,pc_y+d)
    Yaux <- Y[Indexrow,Indexcol]
    tx <- (Indexrow-1)/2^(niveau)
    ty <- (Indexcol-1)/2^(niveau)
    #compute the associated X 
    Xaux <- rbind(as.vector(matrix(rep(tx,length(ty)),nrow=length(ty),byrow=TRUE)),rep(ty,length(tx)),as.vector(t(Yaux)))

if (sum(Xaux[3,]!=0)==0){ 
###################################
       #Simulation for the midle point 
       ###################################
       x<-sqrt(R(cbind(tc,tc)))*rnorm(1)
       X<-cbind(X,c(tc,x))
       Y[pc_x,pc_y]<-x
       ##############################################
       #Simulation for the point at right
       ##############################################
       tr<-tc+c(2^(-niveau),0)
       x<-sqrt(R(cbind(tr,tr)))*rnorm(1)
       X<-cbind(X,c(tr,x))
       Y[(pc_x+1),pc_y]<-x
       ##############################################
       #Simulation for the point above 
       ##############################################
       tr<-tc+c(0,2^(-niveau))
       x<-sqrt(R(cbind(tr,tr)))*rnorm(1)
       X<-cbind(X,c(tr,x))
       Y[pc_x,(pc_y+1)]<-x
       ##############################################
       #Simulation for the point below
       #only for m=1
       ##############################################
       if (m==1){
            tr<-tc+c(0,-2^(-niveau))
            x<-sqrt(R(cbind(tr,tr)))*rnorm(1)
            X<-cbind(X,c(tr,x))
            Y[pc_x,(pc_y-1)]<-x}
       ##############################################
       #Simulation for the point at left 
       #only for l=1
       ##############################################
       if (l==1){
        tr<-tc+c(-2^(-niveau),0)
        x<-sqrt(R(cbind(tr,tr)))*rnorm(1)
        X<-cbind(X,c(tr,x))
        Y[(pc_x-1),pc_y]<-x}


}
else{ 

if (sum(Xaux[3,]!=0)!=1){

Xaux <- Xaux[,Xaux[3,]!=0] 
    
    #////////compute the indexes of neighbors //////////
    normaux<-Xaux[1:2,]-rbind(tc[1]*rep(1,length(Xaux[1,])),tc[2]*rep(1,length(Xaux[1,])))
    normaux<-normaux[1,]^2+normaux[2,]^2
    I <- sort(normaux,index=TRUE)$ix
    I<-I[1:min(nbNeighbor,length(Xaux[1,]))]
    X <- Xaux[,I]
    
    #////////compute filter and variance //////
    x<-cbind(X[1:2,],tc)
    out<-.Call("vf",as.double(x),body(R),new.env())
    var<-out[1]; filtre<-out[2:length(out)];  

    #////////simulation of x(tc) /////////////////////
    x<-sqrt(var)*rnorm(1) - sum(filtre[1:dim(X)[2]]*X[3,]) 
    X<-cbind(X,c(tc,x))
    
    Y[pc_x,pc_y]<-x

    ##############################################
    #Simulation for the point at right
    ##############################################
    tr<-tc+c(2^(-niveau),0)
    
    normaux<-Xaux[1:2,]-rbind(tr[1]*rep(1,length(Xaux[1,])),tr[2]*rep(1,length(Xaux[1,])))
    normaux<-normaux[1,]^2+normaux[2,]^2
    I <- sort(normaux,index=TRUE)$ix
    I<-I[1:min(nbNeighbor,length(Xaux[1,]))]
    X <- Xaux[,I]

    
    #////////compute filter and variance //////
    x<-cbind(X[1:2,],tr)
    out<-.Call("vf",as.double(x),body(R),new.env())
   
    var<-out[1]; filtre<-out[2:length(out)];  
    
    #////////simulation of x(tr) /////////////////////
    x<-sqrt(var)*rnorm(1) - sum(filtre[1:dim(X)[2]]*X[3,]) 
    X<-cbind(X,c(tr,x))
    
    Y[(pc_x+1),pc_y]<-x

    ##############################################
    #Simulation for the point above 
    ##############################################
    tr<-tc+c(0,2^(-niveau))

    normaux<-Xaux[1:2,]-rbind(tr[1]*rep(1,length(Xaux[1,])),tr[2]*rep(1,length(Xaux[1,])))
    normaux<-normaux[1,]^2+normaux[2,]^2
    I <- sort(normaux,index=TRUE)$ix
    I<-I[1:min(nbNeighbor,length(Xaux[1,]))]
    X <- Xaux[,I]

        
    #////////compute filter and variance //////
    x<-cbind(X[1:2,],tr)
    out<-.Call("vf",as.double(x),body(R),new.env())
    var<-out[1]; filtre<-out[2:length(out)];  
    
    #////////simulation of x(tr) /////////////////////
    x<-sqrt(var)*rnorm(1) - sum(filtre[1:dim(X)[2]]*X[3,]) 
    X<-cbind(X,c(tr,x))
    
    Y[pc_x,(pc_y+1)]<-x
    
    
    ##############################################
    #Simulation for the point below
    #only for m=1
    ##############################################
    if (m==1){
    
    tr<-tc+c(0,-2^(-niveau))
    
    normaux<-Xaux[1:2,]-rbind(tr[1]*rep(1,length(Xaux[1,])),tr[2]*rep(1,length(Xaux[1,])))
    normaux<-normaux[1,]^2+normaux[2,]^2
    I <- sort(normaux,index=TRUE)$ix
    I<-I[1:min(nbNeighbor,length(Xaux[1,]))]
    X <- Xaux[,I]

    
    #////////compute filter and variance //////
    x<-cbind(X[1:2,],tr)
    out<-.Call("vf",as.double(x),body(R),new.env())
    var<-out[1]; filtre<-out[2:length(out)];  
    


    #////////simulation of x(tr) /////////////////////
    x<-sqrt(var)*rnorm(1) - sum(filtre[1:dim(X)[2]]*X[3,]) 
    X<-cbind(X,c(tr,x))
    
    Y[pc_x,(pc_y-1)]<-x

    }

    ##############################################
    #Simulation for the point at left 
    #only for l=1
    ##############################################
    if (l==1){
    tr<-tc+c(-2^(-niveau),0)
    
    normaux<-Xaux[1:2,]-rbind(tr[1]*rep(1,length(Xaux[1,])),tr[2]*rep(1,length(Xaux[1,])))
    normaux<-normaux[1,]^2+normaux[2,]^2
    I <- sort(normaux,index=TRUE)$ix
    I<-I[1:min(nbNeighbor,length(Xaux[1,]))]
    X <- Xaux[,I]

    #////////compute filter and variance //////
    x<-cbind(X[1:2,],tr)
    out<-.Call("vf",as.double(x),body(R),new.env())
    var<-out[1]; filtre<-out[2:length(out)];  
    
    #////////simulation of x(tr) /////////////////////
    x<-sqrt(var)*rnorm(1) - sum(filtre[1:dim(X)[2]]*X[3,]) 
    
    Y[(pc_x-1),pc_y]<-x

    }
}
else{
###################################
         #Simulation for the midle point 
         ###################################
         Xaux <- Xaux[,Xaux[3,]!=0]
         tv <- Xaux[1:2]
         alpha <- R(cbind(tc,tv))/R(cbind(tv,tv))
         std <- sqrt(R(cbind(tc,tc))-alpha*R(cbind(tc,tv)))
         x <- std*rnorm(1)+alpha*Xaux[3]
         X<-cbind(X,c(tc,x))
         Y[pc_x,pc_y]<-x
         ##############################################
         #Simulation for the point at right
         ##############################################
         tr<-tc+c(2^(-niveau),0)
         alpha <- R(cbind(tr,tv))/R(cbind(tv,tv))
         std <- sqrt(R(cbind(tr,tr))-alpha*R(cbind(tr,tv)))
         x <- std*rnorm(1)+alpha*Xaux[3]       
         X<-cbind(X,c(tr,x))
         Y[(pc_x+1),pc_y]<-x
         ##############################################
         #Simulation for the point above 
         ##############################################
         tr<-tc+c(0,2^(-niveau))
         alpha <- R(cbind(tr,tv))/R(cbind(tv,tv))
         std <- sqrt(R(cbind(tr,tr))-alpha*R(cbind(tr,tv)))
         x <- std*rnorm(1)+alpha*Xaux[3]       
         X<-cbind(X,c(tr,x))
         Y[pc_x,(pc_y+1)]<-x
         ##############################################
         #Simulation for the point below
         #only for m=1
         ##############################################
         if (m==1){
            tr<-tc+c(0,-2^(-niveau))
            alpha <- R(cbind(tr,tv))/R(cbind(tv,tv))
            std <- sqrt(R(cbind(tr,tr))-alpha*R(cbind(tr,tv)))
            x <- std*rnorm(1)+alpha*Xaux[3]
            X<-cbind(X,c(tr,x))
            Y[pc_x,(pc_y-1)]<-x}
         ##############################################
         #Simulation for the point at left 
         #only for l=1
         ##############################################
         if (l==1){
             tr<-tc+c(-2^(-niveau),0)
             alpha <- R(cbind(tr,tv))/R(cbind(tv,tv))
             std <- sqrt(R(cbind(tr,tr))-alpha*R(cbind(tr,tv)))
             x <- std*rnorm(1)+alpha*Xaux[3]
             X<-cbind(X,c(tr,x))
             Y[(pc_x-1),pc_y]<-x}
}
}
    
}
}
niveau<-niveau+1
Z<-Y
}

time <- proc.time() - time

return(list(Zrow=0:2^(niveau-1),Zcol=0:2^(niveau-1),Z=Z,time=time[[3]]))
}
