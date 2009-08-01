#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP mkansvec3(SEXP x)
{
SEXP ans;
int i,nx;
nx=length(x);
PROTECT(ans=allocVector(REALSXP,nx));
for(i=0;i<nx; i++){
REAL(ans)[i]=REAL(x)[i];
}
UNPROTECT(1);
return ans;
}

double feval3(SEXP x, SEXP func, SEXP rho){
defineVar(install("x"),mkansvec3(x),rho);
return(REAL(eval(func,rho))[0]);
}

SEXP vf3(SEXP x, SEXP R, SEXP rho)
{
double sigma2f, stdaux, beta, r, s, S, t;
int ts, kk, k, j, i;
SEXP ans, u;

ts=length(x)/3;
double Af[ts-1];


double b[ts-2], a[ts-2];
double Abold[ts-1][ts-1];
double Abnew[ts-1][ts-1]; 
double sigma2bold[ts],sigma2bnew[ts];
double tabR[ts];



PROTECT(ans=allocVector(REALSXP,ts));
PROTECT(u=allocVector(REALSXP,6));

REAL(u)[0]=REAL(x)[0];
REAL(u)[1]=REAL(x)[1];
REAL(u)[2]=REAL(x)[2];
REAL(u)[3]=REAL(x)[0];
REAL(u)[4]=REAL(x)[1];
REAL(u)[5]=REAL(x)[2];

sigma2f=feval3(u,R,rho);
sigma2bold[0]=sigma2f;

for (kk=2; kk<=ts; kk++){ 

         REAL(u)[0]=REAL(x)[3*kk-3];
         REAL(u)[1]=REAL(x)[3*kk-2];
         REAL(u)[2]=REAL(x)[3*kk-1];
         REAL(u)[3]=REAL(x)[3*kk-3];
         REAL(u)[4]=REAL(x)[3*kk-2];
         REAL(u)[5]=REAL(x)[3*kk-1];

         sigma2f=feval3(u,R,rho);
        
         for(i=0;i<kk-1;i++){Af[i]=0;
           for(j=0;j<kk-1;j++){Abnew[i][j]=0;}
         }

         /*#******init*****************/

         sigma2bnew[0]=sigma2f;
         tabR[kk-1]=sigma2f;

         for (k=1; k<kk; k++){

    /*       #*****calcul des filtres Af_k *** */
           
           REAL(u)[3]=REAL(x)[3*(kk-k)-3];
           REAL(u)[4]=REAL(x)[3*(kk-k)-2];
           REAL(u)[5]=REAL(x)[3*(kk-k)-1];

           t=feval3(u,R,rho); /*ok*/
           tabR[kk-1-k]=t;
           s=sigma2bold[k-1];
//printf("%f",s);
        

           /*#//// calcul de betadelta(u,v)*/

           if (k==1){
           stdaux=sqrt(sigma2f*s);
           if (stdaux ==0){beta=0;}
           if (stdaux !=0){beta=t/stdaux;}
           }
           if (k!=1){

                for (j=0 ; j<k-1 ;j++){
                   b[j]=Abold[k-2][j]; 
                }
                stdaux=sqrt(sigma2f*s);
           
           if (stdaux ==0){beta=0;}
           if (stdaux !=0){
                  S=0;
                  for (j=0; j<k-1; j++){
                  S=S+b[j]*tabR[kk-1-k+j+1];
                  }
                  beta=(t+S)/stdaux; 
                  }
           }

           if (abs(beta)>1){/*stop("message")*/} 
           
           if (stdaux ==0){
            Af[k-1]=0;
            Abnew[k-1][k-1]=0;}

           if (stdaux !=0){
            r=sqrt(sigma2f/s);
            Af[k-1]=-beta*r;
            Abnew[k-1][k-1]=-beta/r ;}

           if (k>1){
                   
                 for (j=0 ;j<k-1 ;j++){a[j]=Af[j];}
                 for (j=0; j<k-1;j++){
                 Af[j]=a[j]+Af[k-1]*b[k-2-j];
                 Abnew[k-1][j]=b[j]+Abnew[k-1][k-1]*a[k-2-j];}
           }

           sigma2f=(1-pow(beta,2))*sigma2f;
           sigma2bnew[k]=(1-pow(beta,2))*s;
        

           /*#********fin calcul filtre Af ***/
         } /*#end for k

         #stockage*/

         for(j=0;j<ts;j++){
         sigma2bold[j]=sigma2bnew[j];
         for(i=0;i<ts-1;i++){
         Abold[i][j]=Abnew[i][j];
         }
         }
       
         } /*#end for kk*/

//printf("fin");       
REAL(ans)[0]=sigma2f;
for (j=0 ; j<ts-1 ; j++){
REAL(ans)[j+1]=Af[ts-2-j];
}
  
 UNPROTECT(2);
 return(mkansvec3(ans));
}
