printhyper<-function(res){
n <- dim(res$W1)[1]

M=3
x<-seq(-M,M,length=n);
y<-seq(-M,M,length=n);
z<-matrix(0,n,n)
for (i in 1:n){
for (j in 1:n){
z[i,j]<-sqrt(1+(x[i]^2+y[j]^2))
}
}

zlim <- range(cbind(res$W1))
zlen <- zlim[2] - zlim[1] + 1
colramp <- colorRampPalette(brewer.pal(9, "Blues"))

colaux <- colramp(20*zlen)
col1<-colaux[20*(res$W1-zlim[1]+1)]

open3d()
rgl.surface(x,y,z,col=col1,lit=FALSE)
}
