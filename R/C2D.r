C2D<-function(H){(pi^(3/2)*gamma(H+1/2)/(H*sin(pi*H)*gamma(2*H)*gamma(H+1)))^(1/2)}



#C2D<-function(N){
#H<-seq(from=0,to=1,length=2^N+1)

#L<--100
#M<-301
#dx<-2*(-L)/M
#dy<-dx
#Const<-H

#for (i in 2:(2^N+1)){
#S<-0
#G<-H[i]

#for (j in 1:(M-1)){
#for (k in 1:(M-1)){
#S<-S+sin((L+j*dx)/2)^2/((L+j*dx)^2+(L+k*dy)^2)^(1+G)*dx*dy
#}
#}

#Const[i]<-4*S
#}

#Const[0]<-0
#return(Const)
#}
