### H.test.R  (2006-10)
###
###    Performs Hurst parameter test for random gaussian field
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

H.test<-function(Z,H,alternative = "two.sided",conf.level = 0.95){

##    INPUT VARIABLES
#########################
##  Z   : matrix of size NxN 
##        matrix associated with the sample path of one  
##        fractal Brownian field
##        (N=2^(nblevel)+1 and Z[i,j] gives the value of 
##        the process at point (2^(-i*nblevel),2^(-j*nblevel)))  
##  H   : Hurst parameter to be tested
##  alternative = c("two.sided", "less", "greater")
##  conf.level : real between 0 and 1. It represents the confidence level of the interval


##    OUTPUT VARIABLES
##########################
## statistic: the value of the t-statistic.
## p.value: the p-value for the test.
## conf.int: a confidence interval for H
## gamma2 : the constant of the variance used in the test


##  TEST ON INPUT VARIABLES
##############################

if ((is.numeric(Z)==FALSE)||(is.matrix(Z)==FALSE)){
 stop("Message from H.test.R: Z is not of valid type")}

if (dim(Z)[1]!=dim(Z)[2])
 stop("Message from H.test.R: Z is not of valid type")

if ((is.numeric(H)==FALSE)||(is.vector(H)==FALSE))
 stop("Message from H.test.R: H is not of valid type")
 
if ((length(H)!=1)||(H<=0)||(H>1))
 stop("Message from H.test.R: H is not of valid type")
 
if ((alternative!="two.sided")&(alternative!="less")&(alternative!="greater"))
 stop("Message from H.test.R: alternative is not of valid type")

if ((conf.level<=0)||(conf.level>=1))
 stop("Message from H.test.R: conf.level is not of valid type")

gamma2 <- gamma2(H)
hatH <- quadvar(Z)
nblevel <- log((dim(Z)[1]-1),2) 
statistic <- (hatH-H)/sqrt(gamma2)*(2^(nblevel))

if (alternative == "two.sided")       
    p.value <- 2*pnorm(-abs(statistic))
if (alternative == "less") 
    p.value <- pnorm(statistic)
if (alternative == "greater") 
    p.value <- pnorm(-statistic)
    
conf.int <- c(hatH-qnorm(conf.level)*sqrt(gamma2)/(2^(nblevel)),hatH+qnorm(conf.level)*sqrt(gamma2)/(2^(nblevel))) 
return(list(statistic=statistic,p.value=p.value,conf.int=conf.int,gamma2=gamma2))
}
