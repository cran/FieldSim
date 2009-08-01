printsphere<-function(res){

n <- dim(res$W1)[1]
eps <- 1/4+0.01

x<-seq(-1+eps,1-eps,length=n);
y<-seq(-1+eps,1-eps,length=n);
z<-matrix(0,n,n)

for (i in 1:n){
 for (j in 1:n){
   z[i,j]<-sqrt(1-(x[i]^2+y[j]^2))
 }
}

zlim <- range(cbind(res$W1,res$W2,res$W3,res$W4,res$W5,res$W6))
zlen <- zlim[2] - zlim[1] + 1
colramp <- colorRampPalette(brewer.pal(9, "Blues"))

colaux <- colramp(20*zlen)
col1<-colaux[20*(res$W1-zlim[1]+1)]
col2<-colaux[20*(res$W2-zlim[1]+1)]
col3<-colaux[20*(res$W3-zlim[1]+1)]
col4<-colaux[20*(res$W4-zlim[1]+1)]
col5<-colaux[20*(res$W5-zlim[1]+1)]
col6<-colaux[20*(res$W6-zlim[1]+1)]

open3d()
rgl.surface(x,y,z,col=col1,lit=FALSE)
rgl.surface(x,y,-z,col=col2,lit=FALSE)
rgl.surface(x,y,z,coords=c(1,3,2),col=col3,lit=FALSE)
rgl.surface(x,y,-z,coords=c(1,3,2),col=col4,lit=FALSE)
rgl.surface(x,y,z,coords=c(2,1,3),col=col5,lit=FALSE)
rgl.surface(x,y,-z,coords=c(2,1,3),col=col6,lit=FALSE)
}
