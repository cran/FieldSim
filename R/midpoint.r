### midpoint.R  (2006-10)
###
###    fractional Brownian field simulation by the midpoint displacement method
###
### Copyright 2006-10 Alexandre Brouste and Sophie Lambert-Lacroix
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
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

midpoint <- function(H,nblevel=8){

##    INPUT VARIABLES
#########################
##  H   : real in ]0,1[
##      Hurst parameter of the fractal Brownian field
##  nblevel  : positive integer
##      Level associated with the regular space discretization 
##      Form of reguar space discretization : [[0:2^(nblevel)]/2^(nblevel)]^2

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

if ((is.numeric(H)==FALSE)||(H>=1)||(H<=0)){
 stop("Message from midpoint.R: H is not of valid type")}
 
if ((is.numeric(nblevel)==FALSE)||(round(nblevel)-nblevel!=0)||(nblevel<0)){
 stop("Message from midpoint.R: nblevel is not of valid type")}

time <- proc.time()
Z<-matrix(0,2,2)
Z[1,2] <- rnorm(1)
Z[2,1] <- rnorm(1)
Z[2,2] <- rnorm(1)*sqrt(2^H)

level<-1
while (level<=nblevel){
Y<- matrix(0,2^(level)+1,2^(level)+1)


#--------Build Y----------------------

for (l in 1:(2^(level)+1)){ #les colonnes
for (m in 1:(2^(level)+1)){ #les lignes

if (((m/2-floor(m/2))!=0) & ((l/2-floor(l/2))!=0)){
Y[m,l]<-Z[((m-1)/2+1),((l-1)/2+1)] 
}

}}

#--------------Simulation of the other points --------------

    #////////variance of the center point  ////// 
    varC<-(1-1/4*2^(H)-1/8*2^(2*H))*2^(-2*level*H+H)
    
    #////////variance of the remaining points ////// 
    varA<-(1-2^(2*H-2))/2^(2*level*H)
    

for (m in 1:2^(level-1)){ 
for (l in 1:2^(level-1)){ 

    pc_x<-2*l
    pc_y<-2*m

   
    ###################################
    #Simulation of the center point
    ###################################

    Y[pc_x,pc_y]<-sqrt(varC)*rnorm(1) + 1/4*(Y[pc_x-1,pc_y-1]+Y[pc_x-1,pc_y+1]+Y[pc_x+1,pc_y+1]+Y[pc_x+1,pc_y-1]) 
    
    ##############################################
    #Simulation for the point at right
    ##############################################
    
    Y[pc_x+1,pc_y]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x+1,pc_y-1]+Y[pc_x+1,pc_y+1])
 
    ##############################################
    #Simulation for the point above 
    ##############################################
    
    Y[pc_x,pc_y+1]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y+1]+Y[pc_x+1,pc_y+1])
  
    ##############################################
    #Simulation for the point below
    #only for m=1
    ##############################################
    if (m==1){
        
    Y[pc_x,pc_y-1]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y-1]+Y[pc_x+1,pc_y-1])

    }

    ##############################################
    #Simulation for the point at left 
    #only for l=1
    ##############################################
    if (l==1){
    
    Y[pc_x-1,pc_y]<-sqrt(varA)*rnorm(1) + 1/2*(Y[pc_x-1,pc_y-1]+Y[pc_x-1,pc_y+1])
 
    }
   
}
}

level<-level+1
Z<-Y
}

time <- proc.time() - time

return(list(Zrow=0:2^(level-1),Zcol=0:2^(level-1),Z=Z,time=time[[3]]))
}
