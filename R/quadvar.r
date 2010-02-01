### quadvar.R  (2006-10)
###
###    Estimation of the Hurst parameter of the fractal Brownian field
###             by the quadratic variations method
###
### Copyright 2006-10 Alexandre Brouste and Sophie Lambert-Lacroix
###
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

quadvar <- function(manifold,res){

##    INPUT VARIABLES
#########################
##  manifold  : manifold
##  res       : vector of length N constitued of the simulated field


##    OUTPUT VARIABLES
##########################
## H: real  Estimation of the Hurst parameter of the fractal Brownian field


##  TEST ON INPUT VARIABLES
##############################

if(missing(manifold)){ 		
	cat("Error from quadvar.R: parameter manifold is missing\n")
	return(NULL)
}	
	
if(!isS4(manifold)){ 
	cat("Error from quadvar.R: parameter manifold is not of type manifold\n")
	return(NULL)
}else if(!class(manifold)[1]=="manifold"){
	cat("Error from quadvar.R: parameter manifold is not of type manifold\n")
	return(NULL)
}		
	

if(missing(manifold)){ 		
	cat("Error from quadvar.R: parameter res is missing\n")
	return(NULL)
}		
		
if (!is.vector(res)){
	cat("Error from quadvar.R: res is not of valid type")
	return(NULL)
}

	
nameofgrid<-whichgrid(manifold)	

if (manifold@name=="plane"){	
	
	if (nameofgrid=="regular"){
		
		N<-sqrt(length(res))
		Z<-as.matrix(res,N,N)
		
		
	}else if(nameofgrid=="visualization"){
		

		N<-sqrt(length(res))	
		Ng<-log(N-1)/log(2) 	
	

		f<-t(res)	
		
		# Construct the Z matrix
		#########################
		Z<-matrix(f[1:4],2,2,byrow=TRUE)
		indice <- 5
		niveau <- 1
		while (niveau<=Ng){
			Y <- matrix(0,2^(niveau)+1,2^(niveau)+1)
			for (l in 1:(2^(niveau)+1)){ #columns
				for (m in 1:(2^(niveau)+1)){ #rows
					if (((m/2-floor(m/2))!=0) & ((l/2-floor(l/2))!=0))
					Y[m,l]<-Z[((m-1)/2+1),((l-1)/2+1)] 
				}}
			
			for (m in 1:2^(niveau-1)){ 
				for (l in 1:2^(niveau-1)) {
					pc_x<-2*l
					pc_y<-2*m
					Y[pc_x,pc_y]<-f[indice]
					indice<-indice+1
					Y[(pc_x+1),pc_y]<-f[indice]
					indice<-indice+1
					Y[pc_x,(pc_y+1)]<-f[indice] 
					indice<-indice+1
					if (m==1){Y[pc_x,(pc_y-1)]<-f[indice]
						indice<-indice+1}
					if (l==1){Y[(pc_x-1),pc_y]<-f[indice]
						indice<-indice+1}
				}
			}
			niveau<-niveau+1
			Z<-Y
		}		
		
		

	}else{
		cat("Error from quadvar.R: no estimator have been implemented for this grid")	
		return(NULL)
	}

	
	
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
	return(H)
	
}else{
	cat("Error from quadvar.R: no estimator have been implemented for this manifold")	
	return(NULL)	
}
	

}
