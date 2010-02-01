#################################################################
#######                    Fieldsim                      ########
#################################################################

## fieldsim.R  (2006-15)
##
##    Random field simulation by the method FieldSim method
##
## Copyright 2006-15 Alexandre Brouste and Sophie Lambert-Lacroix


##    INPUT VARIABLES
#################################################################
##	manifold	: manifold (R manifold object)
##  R			: covariance function (R function)
##  Ne          : number of points to be simulated with the exact
##                method
##  nbNeighbor  : number of neighbors to be considered in the 
##                fieldsim method
#################################################################



##    OUTPUT VARIABLES
#################################################################
## function returns the value of the process of covariance R on
## the atlas of the manifold with the FieldSim method
#################################################################

fieldsim <-function(manifold,R,Ne,nbNeighbor=4){
	
#-----------------Tests on input variables-----------------------
		
if(missing(manifold)){ 		
	cat("Error from fieldsim.R: parameter manifold is missing\n")
	return(NULL)
}	
	
if(!isS4(manifold)){ 
	cat("Error from fieldsim.R: parameter manifold is not of type manifold\n")
	return(NULL)
}else if(!class(manifold)[1]=="manifold"){
	cat("Error from fieldsim.R: parameter manifold is not of type manifold\n")
	return(NULL)
}
	
mesh<-manifold@atlas	
	

if (!is.function(R)){
	cat("Error from fieldsim.R: parameter R is not of function type\n")
	return(NULL)
}


dimen<-dim(mesh)[1]
		
#Traitement des NaN
		
meshnan<-is.nan(mesh)
indexNan<-meshnan[1,]
for (i in 1:dimen){ 
	indexNan<-indexNan|meshnan[i,]
}
		
meshtmp<-mesh
mesh<-matrix(mesh[1:dimen,!indexNan],nrow=dimen)

Ntmp<-dim(meshtmp)[2]
N<-dim(mesh)[2]	
	

	
if(missing(Ne)){
		Ne<-N
}

if(!is.numeric(Ne)){
	cat("Error from fieldsim.R: Ne must be of numeric type\n")
	return(NULL)

}
	
if (Ne<=0|Ne>N){
	cat("Error from fieldsim.R: Ne have must be up to 0 and less than the size of the atlas\n")
	return(NULL)
}	

if(!is.numeric(nbNeighbor)){
	cat("Error from fieldsim.R: nbNeighbor must be of numeric type\n")
	return(NULL)
		
}
	
if(Ne<=nbNeighbor){
	cat("Error from fieldsim.R: Ne must be bigger than nbNeighbor\n")
	return(NULL)
}	
	

	
	matcov<-matrix(0,Ne,Ne)
		
	for(i in 1:Ne){
		for (j in i:Ne){
			matcov[i, j]<- R(mesh[,i],mesh[,j])		
			matcov[j, i]<- matcov[i, j]
		}
	}

	iso<-apply(matcov==0,2,all)
	whereo<-which(iso)	

	if (any(iso)){
		matcovbis<-matcov[-whereo,-whereo]
		Nbis<-dim(matcovbis)[2]
	}else{
		matcovbis<-matcov
		Nbis<-Ne
	}
			
	L <-matrix(0,Nbis,Nbis)	
	L <- chol(matcovbis)
	Z <- rnorm(Nbis)
	restmp <- t(L) %*% Z
	
	if(any(iso)){
		res<-rep(0,length=Ne)
		res[-whereo]<-as.vector(restmp)
	}else{
		res<-as.vector(restmp)
	}	
	
	
	if (Ne<N){
		
		d<-manifold@distance
		dmat=nbNeighbor+1
				
		for (ts in (Ne+1):N){
			D<-0
			for (k in 1:(ts-1)){
				D[k]<- d(mesh[,ts],mesh[,k])
			}
			
			Voisins<-sort(D,index=TRUE)$ix[1:nbNeighbor]
			xx<-matrix(mesh[1:dimen, c(Voisins,ts)],nrow=dimen)	
					
			Rmat = matrix(0,dmat,dmat)

			for (i in 1:dmat){
				for (j in 1:dmat){ 
						Rmat[i,j]<- R(xx[,i],xx[,j])
				}
			}
			
			ismato<-apply(Rmat==0,2,all)
			wheremato<-which(ismato)	
			
			if (any(ismato)){
				Rmatbis<-Rmat[-wheremato,-wheremato]
			}else{
				Rmatbis<-Rmat
			}
			
			
			Li=t(chol(Rmatbis))  
			D=diag(diag(Li))
			Ld<-dim(D)[2]	
			Lin=Li%*%solve(D)
			Linmoins1 = solve(Lin)
					
			xvec <- res[Voisins]
			if (any(ismato)){
				xsim<-0
			}else{
				xsim<-D[dmat,dmat]*rnorm(1)-sum(Linmoins1[dmat,1:(dmat-1)]*xvec)
			}
				
			res=c(res,xsim)	
			}
		}else{}
	
	
	
	restmp<-rep(NaN,Ntmp)
	restmp[!indexNan]<-res
	res<-restmp	
	return(res)
}


#################################################################
#######                        Fin                         ######
#################################################################