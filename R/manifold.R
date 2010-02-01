setMethod("initialize", "manifold",
function(.Object, name=NULL, atlas=NULL, distance=NULL, origin=NULL){  
	
	if(!is.null(name)){
		.Object@name <- name
	}
	
	if(!is.null(atlas)){ 
		.Object@atlas <- atlas
	}
	
	if(!is.null(distance)){
		.Object@distance <- distance
	}
	
	if(!is.null(origin)){
		.Object@origin <- origin
	}
	
	return(.Object)
})


#################################################################
#######                    Fieldsim                      ########
#################################################################

## setManifold.R  (2006-15)
##
##    
##
## Copyright 2006-15 Alexandre Brouste et Sophie Lambert-Lacroix


##    INPUT VARIABLES
#################################################################
##  name     : name of the manifold (type character)
##  atlas    : atlas of the manifold (type matrix)
##  distance : distance on the manifold (type function)
##  origin	 : origin of the manifold (type matrix)
#################################################################



##    OUTPUT VARIABLES
#################################################################
##   function returns an object manifold  
#################################################################



setManifold <-
function(name, atlas, distance, origin){
	
if(missing(name)){ 		
	cat("Error from setManifold.R: user must name the manifold\n")
	return(NULL)
}
	
if(!is.character(name)){
	cat("Error from setManifold.R: name must be a character string\n")
	return(NULL)
}else{	
#The "plane" manifold
	
	if(name=="plane"){
		
		Ng<-3
		mesh<-NULL
		for (l in 0:1){
			for (m in 0:1){
				mesh <- cbind(mesh,rbind(l,m))  #Grille grossière
			}
		}
		
		niveau <- 1
		while (niveau<=Ng){ #parametre qui donne le nombre de rafinement a faire
			for (m in 1:2^(niveau-1)){ 
				for (l in 1:2^(niveau-1)) {
					pc_x<-2*l
					pc_y<-2*m
					tr_x<-(pc_x-1)/2^(niveau)
					tr_y<-(pc_y-1)/2^(niveau)
					mesh <- cbind(mesh,rbind(tr_x,tr_y))
					mesh <- cbind(mesh,rbind(tr_x+2^(-niveau),tr_y))
					mesh <- cbind(mesh,rbind(tr_x,tr_y+2^(-niveau)))
					if (m==1){mesh <- cbind(mesh,rbind(tr_x,tr_y-2^(-niveau)))}
					if (l==1){mesh <- cbind(mesh,rbind(tr_x-2^(-niveau),tr_y))}
				}
			}
			niveau<-niveau+1
		}
		dimnames(mesh)<-NULL
		
		return(new("manifold",
			   name=name,
			   atlas=as.matrix(mesh),
			   distance=function(xi,xj){return(sqrt(t(xi-xj)%*%(xi-xj)))},
			   origin=rbind(0,0))
		)

	}
	
#The "sphere" manifold	
	
	if(name=="sphere"){
		
		eps<-1/4
		N<-12
		x<-seq(-1+eps,1-eps,length=N);
		z<-matrix(0,N,N)
		
		for (i in 1:N){
			for (j in 1:N){
				suppressWarnings(z[i,j]<-sqrt(1-(x[i]^2+x[j]^2)))
			}
		}
		
		W1<-rbind(rep(x,each=N),rep(x,N),as.vector(z))
		W2<-rbind(rep(x,each=N),rep(x,N),-as.vector(z))
		W3<-rbind(rep(x,each=N),as.vector(z),rep(x,N))
		W4<-rbind(rep(x,each=N),-as.vector(z),rep(x,N))
		W5<-rbind(as.vector(z),rep(x,N),rep(x,each=N))
		W6<-rbind(-as.vector(z),rep(x,N),rep(x,each=N))
		
		atlassphere<-cbind(W1,W2,W3,W4,W5,W6)
		
		return(new("manifold",
				name=name,
				atlas=atlassphere, 
				distance=function(xi,xj){ #Distance on the sphere
				   u <- sum(xi*xj)
				   if (u<(-1))
				   u<--1
				   if (u>1)
				   u<-1
				   return(acos(u))},
				origin=rbind(1,0,0))
		)

	}
	
#The "hyperboloid" manifold
		
	if(name=="hyperboloid"){
		
		M=3
		res<-0
		N<-12
		x<-seq(-M,M,length=N);
		z<-matrix(0,N,N)
		for (i in 1:N){
			for (j in 1:N){
				z[i,j]<-sqrt(1+(x[i]^2+x[j]^2))
			}
		}	
		atlashyper<-rbind(rep(x,each=N),rep(x,N),as.vector(z))
		
		return(new("manifold",
				name=name,
				atlas=atlashyper,
				distance=function(xi,xj){    #Distance on the hyperboloid
				   u <- -xi[1]*xj[1]-xi[2]*xj[2]+xi[3]*xj[3]
				   if (u<1){u<-1}
				   return(acosh(u))},
				origin=rbind(0,0,1))
		)
   }
   
#Users manifolds
	
		
		
		if(missing(atlas)){
			cat("Error from setManifold.R: atlas must be set\n")
			return(NULL)
		}
				   
		if(missing(distance)){
			cat("Error from setManifold.R: distance must be set\n")
			return(NULL)	   
		}
				   
		if(missing(origin)){
			stop("Error from setManifold.R: no origin have been set\n")
		    return(NULL)
		}
		
		if(!is.matrix(atlas)){
			cat("Error from setManifold.R: atlas must be of matrix type\n")
			return(NULL)
		}
	
		if(!is.function(distance)){
			cat("Error from setManifold.R: distance must be a function\n")
			return(NULL)
		}
	
	
		if(!is.matrix(origin)){
			cat("Error from setManifold.R: origin must be of matrix type\n")
			return(NULL)
		}
	
		if(dim(atlas)[1]!=dim(origin)[1]){
			cat("Error from setManifold.R: atlas and origin have not the same dimension\n")
			return(NULL)
		}

		if(dim(atlas)[1]>dim(atlas)[2]){
			cat("Warning from setManifold.R: dimension is greater than the number of points of the atlas\n")
		}
		
		attr(distance,"source")<-NULL
	
		return(new("manifold", 
				name=name, 
				atlas=atlas, 
				distance=distance, 
				origin=origin)
		)
   
   
}
      
}