### quadvar.R  (2006-10)
###
###    Estimation of the Hurst parameter of the fractal Brownian field
###             by the quadratic variations method
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

quadvar <- function(Z){

##    INPUT VARIABLES
#########################
##  Z   : matrix of size NxN 
##        matrix associated with the sample path of one  
##        fractal Brownian field
##        (N=2^(nblevel)+1 and Z[i,j] gives the value of 
##        the process at point (2^(-i*nblevel),2^(-j*nblevel)))  


##    OUTPUT VARIABLES
##########################
## H: real in ]0,1[
##      Estimation of the Hurst parameter of the fractal Brownian field


##  TEST ON INPUT VARIABLES
##############################

if ((is.numeric(Z)==FALSE)||(is.matrix(Z)==FALSE)){
 stop("Message from quadvar.R: Z is not of valid type")}

if (dim(Z)[1]!=dim(Z)[2]){
 stop("Message from quadvar.R: Z is not of valid type")}


N <- dim(Z)[1]
M <- N-1
V2 <- rep(0,2)

for (m in 1:2){
    K <- floor(M/2^(m-1)) 
    j1 <- seq(from=(2*2^(m-1)+1), to=(2^(m-1)*K+1), by=2^(m-1))
    i2 <- 1:N
    j2 <- seq(from=(2^(m-1)+1), to=(2^(m-1)*(K-1)+1), by=2^(m-1))
    j3 <- seq(from=1, to=(2^(m-1)*(K-2)+1), by=2^(m-1))
    Delta21 <- Z[j1,i2]-2*Z[j2,i2]+Z[j3,i2]
    j4 <- 1:(K-1)
    Delta2 <- Delta21[j4,j1]- 2*Delta21[j4,j2]+Delta21[j4,j3]
    V2[m] <- sum(sum((Delta2)^2))
}

H <- log(V2[2]/V2[1])/(2*log(2))+1

return(H)}
