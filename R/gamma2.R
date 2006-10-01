### gamma2.R  (2006-10)
###
###    Compute variance constant in the H-test
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

gamma2<-function(H,N=10){

##    INPUT VARIABLES
#########################
##  H   : Hurst parameter
##  N : number of iteration  


##    OUTPUT VARIABLES
##########################
## gamma2 : variance to be used in the H-test

##Define filters
a <- c(1,-2,1)

##Compute C1
u0 <- 0 
tildeu0 <- 0  

for (i in 1:3) {
 for (iprim in 1:3) {
  for (j in 1:3) {
   for (jprim in 1:3) {
   cprod=a[i]*a[iprim]*a[j]*a[jprim]
   u0 <- u0 + cprod*((i-iprim)^2+(j-jprim)^2)^H
   tildeu0 <- tildeu0 + cprod*((i-iprim/2-1)^2+(j-jprim/2-1)^2)^H 
   }}}}
C1 <- -u0

##Compute C2

C2 <- 2*u0^2
for (l1 in 0:(N-1)){
 for (l2 in 1:(N-1)){
  v1 <- 0
  for (i in 1:3) {
   for (iprim in 1:3) {
    for (j in 1:3) {
     for (jprim in 1:3) {
        cprod=a[i]*a[iprim]*a[j]*a[jprim]
        v1 <- v1+cprod*((l1+i-iprim)^2+(l2+j-jprim)^2)^H
     }}}}
C2 <- C2+8*v1^2
}}

C3 <- 2*tildeu0^2
for (l1 in 1:(N-1)){
 for (l2 in 0:(N-1)){
  v2 <- 0
  v3 <- 0
  v4 <- 0
  v5 <- 0
  for (i in 1:3) {
   for (iprim in 1:3) {
    for (j in 1:3) {
     for (jprim in 1:3) {
        cprod=a[i]*a[iprim]*a[j]*a[jprim]
        v2 <- v2+cprod*((-l1/2+i-iprim/2-1)^2+(-l2/2+j-jprim/2-1)^2)^H 
        v3 <- v3+cprod*((-l1/2+i-iprim/2-1)^2+(l2/2+j-jprim/2-1)^2)^H 
        v4 <- v4+cprod*((l1/2+i-iprim/2-1)^2+(l2/2+j-jprim/2-1)^2)^H 
        v5 <- v5+cprod*((l1/2+i-iprim/2-1)^2+(-l2/2+j-jprim/2-1)^2)^H 
     }}}}
C3 <- C3+2*(v2^2+v3^2+v4^2+v5^2)
}}


d=2
#gamma2<-(C2*(1+2^(-d))-2*C3/2^(-2*H+d))/(log(2)*C1)^2
gamma2<-(C2*(1+2^(-d))-C3*2^(2*H-1))/(log(2)*C1)^2

return(gamma2)
}
