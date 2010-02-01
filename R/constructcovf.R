#################################################################
#######                    Fieldsim                      ########
#################################################################

## constructcovf.R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste et Sophie Lambert-Lacroix


##    INPUT VARIABLES
#################################################################
##  manifold : a manifold (type manifold)
##  typeproc : the type of covariance (type character)
##  H		 : Hurst parameter for fBm
##  F		 : Hurst function parameter for mBm
#################################################################


##    OUTPUT VARIABLES
#################################################################
##  function returns an covariance function of typeproc process 
#################################################################



constructcovf<-function(manifold,typeproc,H,F){

if(missing(manifold)){ 		
	cat("Error from constructcovf.R: parameter manifold is missing\n")
	return(NULL)
}	

if(!isS4(manifold)){ 
	cat("Error from constructcovf.R: parameter manifold is not of type manifold\n")
	return(NULL)
}else if(!class(manifold)[1]=="manifold"){
	cat("Error from constructcovf.R: parameter manifold is not of type manifold\n")
	return(NULL)
}	
	
	
names=c("fBm","mBm")
	
if(missing(typeproc)){ 		
	cat("Error from constructcovf.R: parameter typeproc is missing\n")
	return(NULL)
}		

if(all(typeproc!=names)){
	cat("Error from constructcovf.R: parameter typeproc does not exist\n")
	return(NULL)
}	
	
if (typeproc == "fBm"){
	
	if (missing(H)){
		cat("Warning from constructcovf.R: parameter H has been set to 0.5\n")	
		H<-0.5
	}	
	
	if (!is.numeric(H)){
		cat("Error from constructcovf.R: parameter H must be numeric\n")	
		return(NULL)
	}
	
	if (H>=1|H<=0){
		cat("Error from constructcovf.R: parameter H must belong to (0,1)\n")	
		return(NULL)
	}
	
	Origine<-manifold@origin
	d<-manifold@distance
	
	if (manifold@name=="sphere"&H>0.5){
		cat("There is no fBm on the sphere for H>0.5")
		return(NULL)
	}
	
	R<-function(xi,xj){
		return(1/2*(d(Origine,xi)^{2*H}+d(Origine,xj)^{2*H}-d(xi,xj)^{2*H}))
	}
	return(R)	
}	
	

if (typeproc=="mBm"){
	
	
	
	if (missing(F)){
		cat("Warning from constructcovf.R: parameter F has been set to constant equal to 0.5\n")	
		F<-function(x){return(0.5)}
	}	
	
	if (!is.function(F)){
		cat("Error from constructcovf.R: parameter F must be a function\n")	
		return(NULL)
	}
	
	
	Origine<-manifold@origin
	d<-manifold@distance
	
	
	R<-function(xi,xj){
		H1<-F(xi)
		H2<-F(xj)
		alpha<-1/2*(H1+H2)
		return(C2D(alpha)^2/(2*C2D(H1)*C2D(H2))*(d(Origine,xi)^(2*alpha)+d(Origine,xj)^(2*alpha)-d(xi,xj)^(2*alpha)))
	}			
	return(R)	
}

}
