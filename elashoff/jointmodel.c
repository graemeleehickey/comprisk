/*******************************************/
/*** Author: Ning Li <nli@PHHP.UFL.EDU>  ***/
/*******************************************/

/*** Nov 2008 ***/
/** joint analysis of longitudinal and competing risks data **/




#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include "cholesky.c"


#define array_size 100000
#define kappa 100

static const double xs[] = {
    2.45340708300901249903e-01,    7.37473728545394358719e-01,
    1.23407621539532300786e+00,    1.73853771211658620678e+00,
    2.25497400208927552311e+00,    2.78880605842813048055e+00,
    3.34785456738321632688e+00,    3.94476404011562521040e+00,
    4.60368244955074427298e+00,    5.38748089001123286199e+00
};

static const double ws[] = {
    4.62243669600610089640e-01,    2.86675505362834129720e-01,
    1.09017206020023320014e-01,    2.48105208874636108814e-02,
    3.24377334223786183217e-03,    2.28338636016353967260e-04,
    7.80255647853206369398e-06,    1.08606937076928169398e-07,
    4.39934099227318055366e-10,    2.22939364553415129254e-13
};






double HAZ(const gsl_matrix *H, const double t);

double CH(const gsl_matrix *H, const double t);

double MulVV(const gsl_vector *Z,const gsl_vector *beta);

void MulV(const gsl_vector *Z,gsl_matrix *ZZ);

void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);

void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);

int inv_matrix(gsl_matrix *x_square);

double Abs(const double a, const double b);

int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

double DiffV(const gsl_vector *veca, const gsl_vector *vecb);

double Min(const double t1, const double t2);

void STAT(gsl_matrix *store,int i,double *mean,double *sd);

int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1,const gsl_matrix *preH2,const gsl_matrix *H2);

void TransM(const gsl_matrix *A, gsl_matrix *B);

int Sbeta(gsl_vector *beta, double *sigma, const gsl_matrix *Y, const int p1a);




int EM(
       gsl_vector *beta,
       gsl_matrix *gamma,
       gsl_vector *vee,
       gsl_matrix *H01,
       gsl_matrix *H02,
       double *sigma,
       gsl_matrix *sig,
       const gsl_matrix *Y,
       const gsl_matrix *C,
       const gsl_vector *M1,
       const int p1a,
       const int maxl
       );


int GetCov(
           gsl_matrix *Cov,
           gsl_vector *beta,
           const gsl_matrix *gamma,
           const gsl_vector *vee,
           const gsl_matrix *H01,
           const gsl_matrix *H02,
           const double sigma,
           gsl_matrix *sig,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const int p1a,
           const int maxl
           );


int GetE(
          gsl_vector *FUNU,
          gsl_vector *FUNUS,
          gsl_matrix *FUNB,
          gsl_matrix *FUNBS,
          gsl_matrix *FUNBU,
          gsl_matrix *FUNE,
          gsl_matrix *FUNUSE,
          gsl_matrix *FUNUE,
          gsl_matrix *FUNW,
          gsl_matrix *FUNWB,
          gsl_matrix *FUNWBS,
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl
          );



double Getloglik(
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl
          );


int Diff(
         const gsl_vector *prebeta,
         const gsl_vector *beta,
         const gsl_matrix *pregamma,
         const gsl_matrix *gamma,
         const gsl_vector *prevee,
         const gsl_vector *vee,
         const gsl_matrix *preH01,
         const gsl_matrix *H01,
         const gsl_matrix *preH02,
         const gsl_matrix *H02,
         const double presigma,
         const double sigma,
         const gsl_matrix *presig,
         const gsl_matrix *sig
         );




int main()
{



    int k,n1,p1,p2,maxl,g=2;
    int p1a;   /* dimension of random effects in Y = 1,2,3 */

/*
    k=140; n1=715; p1=8; p2=5; maxl=6; p1a=2; 


*/

    printf("Enter # of subjects in study: \n");
    scanf("%d", &k);
    printf("# of subjects in study is %d\n",k);

    printf("Enter the total # of observations in Y: \n");
    scanf("%d", &n1);
    printf("The total # of observations in Y is %d\n",n1);

    printf("Enter the dimension of fixed effects (including intercept) in Y: \n");
    scanf("%d", &p1);
    printf("The dimension of fixed effects (including intercept) in Y is %d\n",p1);

    printf("Enter the dimension of fixed effects in C: \n");
    scanf("%d", &p2);
    printf("The dimension of fixed effects in C is %d\n",p2);


    printf("Enter the maximum number of observations per subject in Y: \n");
    scanf("%d", &maxl);
    printf("The maximum number of observations per subject in Y is %d\n",maxl);


    printf("Enter the dimension of random effects in Y: \n");
    scanf("%d", &p1a);
    printf("The dimension of random effects in Y is %d\n",p1a);



    /* allocate space for data */
    gsl_matrix *C = gsl_matrix_alloc(k,p2+2);
    gsl_matrix *Y= gsl_matrix_alloc(n1, p1+p1a+1);    
    gsl_vector *M1= gsl_vector_alloc(k);



    /* read Y matrix  */
    {  
      FILE * f = fopen("y.txt", "r");  
      gsl_matrix_fscanf (f, Y);
      fclose (f);
    } 
   

    /* read C matrix  */
    {  
      FILE * f = fopen("c.txt", "r");
      gsl_matrix_fscanf (f, C);
      fclose (f);
    }


    /* read M1 vector  */
    {  
      FILE * f = fopen("m.txt", "r");
      gsl_vector_fscanf (f, M1);
      fclose (f);
    }


    int i,j,iter,status;


    gsl_matrix * Cov=gsl_matrix_alloc(p1+(p2+1)*g+(p1a+1)*(p1a+2)/2,p1+(p2+1)*g+(p1a+1)*(p1a+2)/2);


    /* allocate space for estimated parameters */    
    gsl_matrix * gamma=gsl_matrix_alloc(g, p2);              
    gsl_vector * vee=gsl_vector_alloc(g-1);                   
    gsl_vector * beta=gsl_vector_alloc(p1);                    

    double sigma;
    gsl_matrix *sig = gsl_matrix_alloc(p1a+1,p1a+1);



    /* allocate space for pre parameters */            
    
    gsl_vector * prebeta=gsl_vector_alloc(p1);                
    gsl_matrix * pregamma=gsl_matrix_alloc(g, p2);                
    gsl_vector * prevee=gsl_vector_alloc(g-1);                                          

    double presigma;
    gsl_matrix *presig = gsl_matrix_alloc(p1a+1,p1a+1);



    /* allocate space for standard error estimates */  

    gsl_vector * vbeta=gsl_vector_alloc(p1);                
    gsl_matrix * vgamma=gsl_matrix_alloc(g, p2);                
    gsl_vector * vvee=gsl_vector_alloc(g-1);   

    double v_sigma;
    gsl_vector *vsig = gsl_vector_alloc((p1a+1)*(p1a+2)/2);


    FILE *output_F; 

    output_F=fopen("output.dat","w");
              if (output_F == NULL)
              {printf("Can't open file\n");
              }





    
    gsl_matrix * FH01 = gsl_matrix_alloc(2,array_size);                 
    gsl_matrix * FH02 = gsl_matrix_alloc(2,array_size);                
    
    gsl_matrix_set_zero(FH01);
    gsl_matrix_set_zero(FH02);


    /* count # events for risk 1 */
    int t, u, v, a=0, b=0;

        for(j=0;j<k;j++)
        {
            u=0;
            if(gsl_matrix_get(C,j,1)==1) 
            {
                if(a>=1)
                {
                    for(t=0;t<a;t++)
                    {
                        if(gsl_matrix_get(FH01,0,t)==gsl_matrix_get(C,j,0))
                        {
                            gsl_matrix_set(FH01,1,t,gsl_matrix_get(FH01,1,t)+1);
                            u=1;
                        }
                    }
                    if(u==0)
                    {
                        t=-1;

                        do
                        {
                            t=t+1;
                        }while(t<=a-1 && gsl_matrix_get(FH01,0,t) < gsl_matrix_get(C,j,0));

                        if(t==a)
                        {
                            gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                            gsl_matrix_set(FH01,1,t,1);
                        }

                        if(t<a)
                        {
                            for(v=a-1;v>=t;v--)
                            {
                                gsl_matrix_set(FH01,0,v+1,gsl_matrix_get(FH01,0,v));
                                gsl_matrix_set(FH01,1,v+1,gsl_matrix_get(FH01,1,v));
                            }
                                gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                                gsl_matrix_set(FH01,1,t,1);
                        }
                      
                        a=a+1;
                    }
                }
               
                if(a==0)
                {
                    gsl_matrix_set(FH01,0,0,gsl_matrix_get(C,j,0));
                    gsl_matrix_set(FH01,1,0,1);
                    a=a+1;
                }
            }
        }


    if(a==0) 
    {
        printf("No failure time information for risk 1; Program exits\n");
        return 0;
    }
     

    gsl_matrix * H01 = gsl_matrix_alloc(3,a);                 /* baseline hazard function for competing risk 1 */

    for(i=0;i<3;i++)
    {
        if(i<=1)
        {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j, gsl_matrix_get(FH01,i,j));
        }
        if(i==2)   
        {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j,0.0001);
        }
    }
                                     

     
    /* count # events for risk 2 */

        for(j=0;j<k;j++)
        {
            u=0;
            if(gsl_matrix_get(C,j,1)==2) 
            {
                if(b>=1)
                {
                    for(t=0;t<b;t++)
                    {
                        if(gsl_matrix_get(FH02,0,t)==gsl_matrix_get(C,j,0))
                        {
                            gsl_matrix_set(FH02,1,t,gsl_matrix_get(FH02,1,t)+1);
                            u=1;
                        }
                    }
                    if(u==0)
                    {
                        t=-1;

                        do
                        {
                            t=t+1;
                        }while(t<=b-1 && gsl_matrix_get(FH02,0,t) < gsl_matrix_get(C,j,0));

                        if(t==b)
                        {
                            gsl_matrix_set(FH02,0,t,gsl_matrix_get(C,j,0));
                            gsl_matrix_set(FH02,1,t,1);
                        }

                        if(t<b)
                        {
                            for(v=b-1;v>=t;v--)
                            {
                                gsl_matrix_set(FH02,0,v+1,gsl_matrix_get(FH02,0,v));
                                gsl_matrix_set(FH02,1,v+1,gsl_matrix_get(FH02,1,v));
                            }
                                gsl_matrix_set(FH02,0,t,gsl_matrix_get(C,j,0));
                                gsl_matrix_set(FH02,1,t,1);
                        }
                        b=b+1;
                    }
                }
               
                if(b==0)
                {
                    gsl_matrix_set(FH02,0,0,gsl_matrix_get(C,j,0));
                    gsl_matrix_set(FH02,1,0,1);
                    b=b+1;
                }
            }
        }
 

    if(b==0) 
    {
        printf("No failure time information for risk 2; Program exits\n");
        return 0;
    }


    gsl_matrix * H02 = gsl_matrix_alloc(3,b);                 /* baseline hazard function for competing risk 2 */

    for(i=0;i<3;i++)
    {
        if(i<=1)
        {
            for(j=0;j<b;j++)    gsl_matrix_set(H02,i,j, gsl_matrix_get(FH02,i,j));
        }
        if(i==2)   
        {
            for(j=0;j<b;j++)    gsl_matrix_set(H02,i,j,0.0001);
        }
    }


    gsl_matrix * preH01 = gsl_matrix_alloc(3,a);               
    gsl_matrix * preH02 = gsl_matrix_alloc(3,b);   




    /* initialize the parameters */

    gsl_matrix_set_zero(gamma);
    gsl_vector_set_zero(vee);
    gsl_vector_set_zero(beta);

    Sbeta(beta,&sigma,Y,p1a);
    gsl_vector_set(vee,0,0);
    gsl_matrix_set_identity(sig); 



    double loglike;


    iter=0;
    do
    {
        iter=iter+1;

        /* store the pre-information */

        gsl_vector_memcpy(prebeta, beta);
        gsl_vector_memcpy(prevee, vee);  
        gsl_matrix_memcpy(pregamma, gamma);
        gsl_matrix_memcpy(preH01, H01);
        gsl_matrix_memcpy(preH02, H02); 
        gsl_matrix_memcpy(presig, sig); 

        presigma=sigma; 
 
        
        /* get new parameter estimates */


        status = EM(beta,gamma,vee,H01,H02,&sigma,sig,Y,C,M1,p1a,maxl);

        printf("iter=%d   status=%d\n",iter,status);

 
            printf("Beta = \n");
            for (i=0;i<p1;i++)
            {
                printf("%f     ", gsl_vector_get(beta,i)); 
                      
            }
            printf("\n");


            printf("Gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                    printf("%f    ", gsl_matrix_get(gamma,i,j));
                }
                printf("\n");
            }
    
            printf("Vee = \n");
            for (i=0;i<g-1;i++)
            {
                printf("%f     ", gsl_vector_get(vee,i)); 
                      
            }
            printf("\n"); 


            printf("sigma = %f\n",sigma);



            printf("Sig = \n");
            for (i=0;i<p1a+1;i++)
            {
                for(j=0;j<p1a+1;j++)
                {
                    printf("%f    ", gsl_matrix_get(sig,i,j));
                }
                printf("\n");
            }

    }while(Diff(prebeta,beta,pregamma,gamma,prevee,vee,preH01,H01,preH02,H02,presigma,sigma,presig,sig)==1
           && status != 100 && iter<100000);

    if(status==100) 
    {
        printf("program stops because of error\n");
        return 0;
    }
    

    if(iter==100000) 
    {
        printf("program stops because of nonconvergence\n");
        return 0;
    }


    if(status != 100 && iter<100000)
    {

        /* if algorithm coverges, compute the variance-covariance matrix of parameters ***/




        status = GetCov(Cov,beta,gamma,vee,H01,H02,sigma,sig,Y,C,M1,p1a,maxl);

        if(status==100) 
        {
            printf("program stops because of error\n");
            return 0;
        }


        if(status != 100)
        {
            status=inv_matrix(Cov);

            if(status==100) 
            {
                printf("program stops because of error\n");
                return 0;
            }

            if(status!=100)
            {

            fprintf(output_F,"Variance-covariance matrix for all the parameters: \n");
            for (i=0;i<Cov->size1;i++)
            {
                for(j=0;j<Cov->size1;j++)
                {
                    fprintf(output_F,"%f    ", gsl_matrix_get(Cov,i,j));
                }
                fprintf(output_F,"\n");
            }


            for (i=0;i<p1;i++)  gsl_vector_set(vbeta,i,gsl_matrix_get(Cov,i,i));

            v_sigma=gsl_matrix_get(Cov,p1,p1);

            for (i=p1+1;i<p1+p2+1;i++)  gsl_matrix_set(vgamma,0,i-p1-1,gsl_matrix_get(Cov,i,i));
            for (i=p1+p2+1;i<p1+2*p2+1;i++)  gsl_matrix_set(vgamma,1,i-p1-p2-1,gsl_matrix_get(Cov,i,i));

            for (i=p1+2*p2+1;i<p1+2*p2+g;i++)  gsl_vector_set(vvee,i-p1-2*p2-1,gsl_matrix_get(Cov,i,i));

            for (i=p1+2*p2+g;i<p1+2*p2+g+(p1a+1)*(p1a+2)/2;i++)  gsl_vector_set(vsig,i-p1-2*p2-g,gsl_matrix_get(Cov,i,i));



            fprintf(output_F,"Beta = \n");
            for (i=0;i<p1;i++)
            {
                fprintf(output_F,"%f     ", gsl_vector_get(beta,i));

            }
            fprintf(output_F,"\n");

            fprintf(output_F,"standard error of beta = \n");
            for (i=0;i<p1;i++)
            {
                fprintf(output_F,"%f     ", sqrt(gsl_vector_get(vbeta,i)));

            }
            fprintf(output_F,"\n");


            fprintf(output_F,"Gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                    fprintf(output_F,"%f    ", gsl_matrix_get(gamma,i,j));
                }  
                fprintf(output_F,"\n");
            }


            fprintf(output_F,"standard error of gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                    fprintf(output_F,"%f    ", sqrt(gsl_matrix_get(vgamma,i,j)));
                }  
                fprintf(output_F,"\n");
            }


            fprintf(output_F,"v = \n");
            for (i=0;i<g-1;i++)
            {
                fprintf(output_F,"%f     ", gsl_vector_get(vee,i));
 
            }
            fprintf(output_F,"\n");


            fprintf(output_F,"standard error of v = \n");
            for (i=0;i<g-1;i++)
            {
                fprintf(output_F,"%f     ", sqrt(gsl_vector_get(vvee,i)));
 
            }
            fprintf(output_F,"\n");



            fprintf(output_F,"sigma2 = %f\n", sigma);
            fprintf(output_F,"standard error of sigma2 = %f\n", sqrt(v_sigma));


            fprintf(output_F,"Sigma = \n");
            for (i=0;i<p1a+1;i++)  
            {
                for (j=0;j<p1a+1;j++)  fprintf(output_F,"%f     ", gsl_matrix_get(sig,i,j));
                fprintf(output_F,"\n");
            }

            fprintf(output_F,"standard error of Sigma = \n");
            for (i=0;i<(p1a+1)*(p1a+2)/2;i++)  fprintf(output_F,"%f\n", sqrt(gsl_vector_get(vsig,i)));


            }
            
        }      
    }

    gsl_matrix_free(FH01);
    gsl_matrix_free(FH02);
    gsl_matrix_free(H01);
    gsl_matrix_free(H02);
    gsl_matrix_free(preH01);
    gsl_matrix_free(preH02);

    gsl_matrix_free(Y);
    gsl_matrix_free(C);
    gsl_vector_free(M1);

    gsl_matrix_free(gamma);
    gsl_vector_free(beta);
    gsl_vector_free(vee);

    gsl_matrix_free(pregamma);
    gsl_vector_free(prebeta);
    gsl_vector_free(prevee);


    gsl_matrix_free(vgamma);
    gsl_vector_free(vbeta);
    gsl_vector_free(vvee);

    gsl_matrix_free(sig);
    gsl_matrix_free(presig);
    gsl_vector_free(vsig);

    gsl_matrix_free(Cov); 

    return 0;


}




int EM(
       gsl_vector *beta,
       gsl_matrix *gamma,
       gsl_vector *vee,
       gsl_matrix *H01,
       gsl_matrix *H02,
       double *sigma,
       gsl_matrix *sig,
       const gsl_matrix *Y,
       const gsl_matrix *C,
       const gsl_vector *M1,
       const int p1a,
       const int maxl
       )

{

    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;    
    
    int n1 = Y->size1;
    int k = M1->size;

    int p,q,j,t,u,i;


    gsl_vector *FUNU=gsl_vector_alloc(k);
    gsl_vector *FUNUS=gsl_vector_alloc(k);
    gsl_matrix *FUNB=gsl_matrix_alloc(p1a,k);  
    gsl_matrix *FUNBS=gsl_matrix_alloc(p1a*(p1a+1)/2,k);  
    gsl_matrix *FUNBU=gsl_matrix_alloc(p1a,k);
   

    gsl_matrix *FUNE=gsl_matrix_alloc(g,k), 
               *FUNUSE=gsl_matrix_alloc(g-1,k),
               *FUNUE=gsl_matrix_alloc(g-1,k);


    gsl_matrix *FUNW=gsl_matrix_alloc(maxl,k), 
               *FUNWB=gsl_matrix_alloc(maxl*p1a,k),
               *FUNWBS=gsl_matrix_alloc(maxl*p1a*(p1a+1)/2,k);


    int status;
    status = GetE(FUNU,FUNUS,FUNB,FUNBS,FUNBU,FUNE,FUNUSE,FUNUE,FUNW,FUNWB,FUNWBS,beta,gamma,vee,H01,H02,
                  *sigma,sig,Y,C,M1,p1a,maxl);



    if (status==100) return status;

    gsl_vector * SX = gsl_vector_alloc(p2);
    gsl_matrix * SXX = gsl_matrix_alloc(p2,p2);
    gsl_matrix * XX = gsl_matrix_alloc(p2,p2);
    gsl_vector * X = gsl_vector_alloc(p2);
    gsl_vector * gammai = gsl_vector_alloc(p2);
  
    double scalef;
  
    gsl_vector * Z = gsl_vector_alloc(p1);                         /* covariates for Y */
    gsl_vector * SZ = gsl_vector_alloc(p1);
    gsl_matrix * SZZ = gsl_matrix_alloc(p1,p1);                    /* store sum of ZZ' */
    gsl_matrix * ZZ = gsl_matrix_alloc(p1,p1);                     /* store ZZ' */
    gsl_vector * xtilde = gsl_vector_alloc(p1a);
    gsl_vector * xtilde1 = gsl_vector_alloc(p1a);

    gsl_vector * bi = gsl_vector_alloc(p1a);
    gsl_matrix * bs = gsl_matrix_alloc(p1a,p1a);



    /* calculate beta, sig */

    gsl_matrix_set_zero(SZZ);
    gsl_vector_set_zero(SZ);

    p=0;
    gsl_matrix_set_zero(sig);

    for(j=0;j<k;j++)
    {
        u=(int)gsl_vector_get(M1,j);

        for(q=p;q<u+p;q++)
        {
            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
            MulV(Z,ZZ);
            gsl_matrix_scale(ZZ,gsl_matrix_get(FUNW,q-p,j));
            gsl_matrix_add(SZZ,ZZ);

            for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNWB,(q-p)*p1a+t,j));
            for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
            gsl_vector_scale(Z,gsl_matrix_get(FUNW,q-p,j)*gsl_matrix_get(Y,q,0)-MulVV(bi,xtilde));
            gsl_vector_add (SZ,Z);
        }
        
        p=p+u;

        gsl_matrix_set(sig,p1a,p1a,gsl_matrix_get(sig,p1a,p1a)+gsl_vector_get(FUNUS,j));
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,t,gsl_matrix_get(sig,t,t)+gsl_matrix_get(FUNBS,t,j));   
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,p1a,gsl_matrix_get(sig,t,p1a)+gsl_matrix_get(FUNBU,t,j));     
        
        for(q=1;q<p1a;q++)
        {
            for(t=0;t<p1a-q;t++)   gsl_matrix_set(sig,t,q+t,gsl_matrix_get(sig,t,q+t)+gsl_matrix_get(FUNBS,p1a+t+(q-1)*(p1a-1),j));
        }

    }

    gsl_matrix_scale(sig,1/(double)k);
    

    status=inv_matrix(SZZ);
    if(status==100)    return status;
    MulM(SZZ,SZ,beta);
  




    /* calculate sigma */

    *sigma=0;
    p=0;

    for(j=0;j<k;j++)
    {
        u=(int)gsl_vector_get(M1,j);

        for(q=p;q<u+p;q++)   
        {
            for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNWB,(q-p)*p1a+t,j));
            for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNWBS,(q-p)*p1a*(p1a+1)/2+t,j));    
        
            for(i=1;i<p1a;i++)
            {
                for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNWBS,(q-p)*p1a*(p1a+1)/2+p1a+t+(i-1)*(p1a-1),j));
            }

            for(t=0;t<p1a;t++)
            {
                for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
            }


            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));   
            for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));

            MulM(bs,xtilde,xtilde1);

            *sigma+=gsl_matrix_get(FUNW,q-p,j)*gsl_pow_2(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))-2*(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))
                    *MulVV(bi,xtilde)+MulVV(xtilde,xtilde1); 


        }
        p=p+u;
    }

    *sigma=*sigma/(double)n1;



    /* calculate H01 H02 */


    double dem, num;

    for(p=0;p<a;p++)
    {
        dem=0;
   
        for(j=0;j<k;j++)
        {
            if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H01,0,p))
            {   
                for(u=0;u<p2;u++)  
                {
                    gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                }

                dem=dem+gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));
                
            }
        }

        gsl_matrix_set(H01,2,p,gsl_matrix_get(H01,1,p)/dem);

    }

     
    for(p=0;p<b;p++)
    {
        dem=0;
   
        for(j=0;j<k;j++)
        {
            if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H02,0,p))
            {   
                for(u=0;u<p2;u++)  
                {
                    gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
                }

                dem=dem+gsl_matrix_get(FUNE,1,j)*exp(MulVV(X,gammai));
                
            }
        }

        gsl_matrix_set(H02,2,p,gsl_matrix_get(H02,1,p)/dem);

    }


   
    /* calculate gamma */

    gsl_matrix_set_zero(SXX);
    gsl_vector_set_zero(SX);
  
  
    for(j=0;j<k;j++)
    {

        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
        }
        
        MulV(X,XX);

        scalef=CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);
        gsl_matrix_scale(XX,scalef); 

        if((int)gsl_matrix_get(C,j,1) == 1)
        {
            scalef=1-CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);  
            gsl_vector_scale(X,scalef); 
        }

        if((int)gsl_matrix_get(C,j,1) != 1)
        {
            scalef=0-CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);  
            gsl_vector_scale(X,scalef); 
        }

        gsl_matrix_add(SXX, XX);
        gsl_vector_add(SX, X);

    }
       

    status=inv_matrix(SXX);
    if(status==100)  return status;
    MulM(SXX,SX,X);
    
    for(j=0;j<p2;j++)   gsl_matrix_set(gamma,0,j,gsl_matrix_get(gamma,0,j)+gsl_vector_get(X,j));





    gsl_matrix_set_zero(SXX);
    gsl_vector_set_zero(SX);
  
    for(j=0;j<k;j++)
    {
        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
        }
        
        MulV(X,XX);

        scalef=CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);
        gsl_matrix_scale(XX,scalef); 

        if((int)gsl_matrix_get(C,j,1) == 2)
        {
            scalef=1-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);  
            gsl_vector_scale(X,scalef); 
        }

        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            scalef=0-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);  
            gsl_vector_scale(X,scalef); 
        }

        gsl_matrix_add(SXX, XX);
        gsl_vector_add(SX, X);

    }
       

    status=inv_matrix(SXX);
    if(status==100)  return status;
    MulM(SXX,SX,X);
    
    for(j=0;j<p2;j++)   gsl_matrix_set(gamma,1,j,gsl_matrix_get(gamma,1,j)+gsl_vector_get(X,j));





    /* calculate vee */

    dem=0; num=0;
    for(j=0;j<k;j++)
    {
        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
        }

        dem=dem+CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUSE,0,j);

        if((int)gsl_matrix_get(C,j,1)==2) 
        {
            num=num+gsl_vector_get(FUNU,j)-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUE,0,j);
        }

        if((int)gsl_matrix_get(C,j,1)!=2) 
        {
            num=num-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUE,0,j);
        }
 
    }

    gsl_vector_set(vee,0,gsl_vector_get(vee,0)+num/dem);



    gsl_vector_free(Z);
    gsl_vector_free(SZ);
    gsl_vector_free(X);
    gsl_vector_free(gammai);
    gsl_vector_free(SX);
    gsl_matrix_free(SZZ);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(SXX);
    gsl_matrix_free(XX);

    gsl_vector_free(bi);
    gsl_vector_free(xtilde);
    gsl_vector_free(xtilde1);
    gsl_matrix_free(bs);


    gsl_vector_free(FUNU);
    gsl_vector_free(FUNUS);
    gsl_matrix_free(FUNB);
    gsl_matrix_free(FUNBS);
    gsl_matrix_free(FUNBU);
         

    gsl_matrix_free(FUNE); 
    gsl_matrix_free(FUNUSE); 
    gsl_matrix_free(FUNUE); 

    gsl_matrix_free(FUNW); 
    gsl_matrix_free(FUNWB); 
    gsl_matrix_free(FUNWBS); 

    return 0;

}



int GetCov(
           gsl_matrix *Cov,
           gsl_vector *beta,
           const gsl_matrix *gamma,
           const gsl_vector *vee,
           const gsl_matrix *H01,
           const gsl_matrix *H02,
           const double sigma,
           gsl_matrix *sig,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const int p1a,
           const int maxl
           )
{

    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;      
    int d =Cov->size1;

    int n1 = Y->size1;
    int k = M1->size;


    int p,q,j,t,u,i,r;


    double temp;

    gsl_vector *S = gsl_vector_alloc(d);
    gsl_matrix *SS= gsl_matrix_alloc(d,d);
    gsl_vector *TS = gsl_vector_alloc(d);


    gsl_vector *FUNU=gsl_vector_alloc(k);
    gsl_vector *FUNUS=gsl_vector_alloc(k);
    gsl_matrix *FUNB=gsl_matrix_alloc(p1a,k);  
    gsl_matrix *FUNBS=gsl_matrix_alloc(p1a*(p1a+1)/2,k);  
    gsl_matrix *FUNBU=gsl_matrix_alloc(p1a,k);
 
   

    gsl_matrix *FUNE=gsl_matrix_alloc(g,k), 
               *FUNUSE=gsl_matrix_alloc(g-1,k),
               *FUNUE=gsl_matrix_alloc(g-1,k);


    gsl_matrix *FUNW=gsl_matrix_alloc(maxl,k), 
               *FUNWB=gsl_matrix_alloc(maxl*p1a,k),
               *FUNWBS=gsl_matrix_alloc(maxl*p1a*(p1a+1)/2,k);


    int status;
    status = GetE(FUNU,FUNUS,FUNB,FUNBS,FUNBU,FUNE,FUNUSE,FUNUE,FUNW,FUNWB,FUNWBS,beta,gamma,vee,H01,H02,
                  sigma,sig,Y,C,M1,p1a,maxl);
    if (status==100) return status;



    gsl_vector * X = gsl_vector_alloc(p2);                         /* covariates for C */
    gsl_vector * RX= gsl_vector_alloc(p2);                         /* covariates for subjects in risk set */
    gsl_vector *SX = gsl_vector_alloc(p2);           
    gsl_vector *SRX= gsl_vector_alloc(p2);
    gsl_vector * gammai = gsl_vector_alloc(p2);

    gsl_vector * Z = gsl_vector_alloc(p1);                         /* covariates for Y */
    gsl_vector * SZ = gsl_vector_alloc(p1);
    gsl_matrix * SZZ = gsl_matrix_alloc(p1,p1);                    /* store sum of ZZ' */
    gsl_matrix * ZZ = gsl_matrix_alloc(p1,p1);                     /* store ZZ' */


    gsl_vector * xtilde = gsl_vector_alloc(p1a);
    gsl_vector * xtilde1 = gsl_vector_alloc(p1a);
    gsl_vector * bi = gsl_vector_alloc(p1a);
    gsl_matrix * bs = gsl_matrix_alloc(p1a,p1a);


    gsl_matrix * VC = gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_matrix * VI = gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_matrix * HE = gsl_matrix_alloc(p1a+1,p1a+1);

    gsl_matrix_memcpy(VC,sig);

    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }
  
    status=inv_matrix(VC);
    if(status==100) return 100;



    /* calculate beta, sig */

    gsl_matrix_set_zero(SZZ);
    gsl_vector_set_zero(SZ);

    p=0;
    gsl_matrix_set_zero(sig);

    for(j=0;j<k;j++)
    {
        u=(int)gsl_vector_get(M1,j);

        for(q=p;q<u+p;q++)
        {
            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
            MulV(Z,ZZ);
            gsl_matrix_scale(ZZ,gsl_matrix_get(FUNW,q-p,j));
            gsl_matrix_add(SZZ,ZZ);

            for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNWB,(q-p)*p1a+t,j));
            for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
            gsl_vector_scale(Z,gsl_matrix_get(FUNW,q-p,j)*gsl_matrix_get(Y,q,0)-MulVV(bi,xtilde));
            gsl_vector_add (SZ,Z);
        }
        
        p=p+u;

        gsl_matrix_set(sig,p1a,p1a,gsl_matrix_get(sig,p1a,p1a)+gsl_vector_get(FUNUS,j));
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,t,gsl_matrix_get(sig,t,t)+gsl_matrix_get(FUNBS,t,j));   
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,p1a,gsl_matrix_get(sig,t,p1a)+gsl_matrix_get(FUNBU,t,j));     
        
        for(q=1;q<p1a;q++)
        {
            for(t=0;t<p1a-q;t++)   gsl_matrix_set(sig,t,q+t,gsl_matrix_get(sig,t,q+t)+gsl_matrix_get(FUNBS,p1a+t+(q-1)*(p1a-1),j));
        }

    }

    gsl_matrix_scale(sig,1/(double)k);
    

    status=inv_matrix(SZZ);
    if(status==100)    return status;
    MulM(SZZ,SZ,beta);
  
  


    p=0;
    gsl_matrix_set_zero(Cov);
    gsl_vector_set_zero(TS);

    for(j=0;j<k;j++)
    {
    
        gsl_vector_set_zero(S);
        u=(int)gsl_vector_get(M1,j);


        /* calculate score for beta */

        gsl_vector_set_zero(SZ); 
   
        for(q=p;q<u+p;q++)
        {
            for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNWB,(q-p)*p1a+t,j));
            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
            for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
            temp = MulVV(beta,Z);

            gsl_vector_scale (Z,gsl_matrix_get(FUNW,q-p,j)*(gsl_matrix_get(Y,q,0)-temp)-MulVV(bi,xtilde));
            gsl_vector_add (SZ,Z);
        }

        gsl_vector_scale (SZ,1/sigma);
        for(q=0;q<p1;q++)  gsl_vector_set(S,q,gsl_vector_get(SZ,q));



        /* calculate score for sigma */

        temp=0;
        for(q=p;q<u+p;q++)   
        {  
            for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNWB,(q-p)*p1a+t,j));
            for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNWBS,(q-p)*p1a*(p1a+1)/2+t,j));    
        
            for(i=1;i<p1a;i++)
            {
                for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNWBS,(q-p)*p1a*(p1a+1)/2+p1a+t+(i-1)*(p1a-1),j));
            }

            for(t=0;t<p1a;t++)
            {
                for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
            }


            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));   
            for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));

            MulM(bs,xtilde,xtilde1);

            temp+=gsl_matrix_get(FUNW,q-p,j)*gsl_pow_2(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))-2*(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))
                    *MulVV(bi,xtilde)+MulVV(xtilde,xtilde1); 

        }

        temp=temp/(2*sigma*sigma);
        temp=temp-(double)u/(2*sigma);
        gsl_vector_set(S,p1,temp);

        p=p+u;




        /* calculate score for gamma */
         

        /*  gamma11, gamma12 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));

        gsl_vector_set_zero(SRX);
        for(r=0;r<a;r++)
        {
            if(gsl_matrix_get(H01,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                gsl_vector_set_zero(SX);
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H01,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                        gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                        gsl_vector_add(SX,RX);                 
                    }
                }

                gsl_vector_scale(SX, gsl_matrix_get(H01,1,r)/(temp*temp));
                gsl_vector_add(SRX,SX);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SX, CH(H01,gsl_matrix_get(C,j,0)));

        gsl_vector_sub(SRX,SX);

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SRX, exp(MulVV(SX,gammai))*gsl_matrix_get(FUNE,0,j));


        if((int)gsl_matrix_get(C,j,1) != 1)
        {
            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(SRX,q));
        }

        if((int)gsl_matrix_get(C,j,1) ==1)
        {
            
            gsl_vector_set_zero(SX);
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                    gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                    gsl_vector_add(SX,RX);
             
                }
            }
            gsl_vector_scale(SX, 1/temp);

            for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));          
            gsl_vector_sub(X,SX); 
             
            gsl_vector_add(X, SRX);

            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(X,q));

        }



        /*  gamma21, gamma22 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));

        gsl_vector_set_zero(SRX);
        for(r=0;r<b;r++)
        {
            if(gsl_matrix_get(H02,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                gsl_vector_set_zero(SX);
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H02,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                        gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q));
                        gsl_vector_add(SX,RX);                 
                    }
                }

                gsl_vector_scale(SX, gsl_matrix_get(H02,1,r)/(temp*temp));
                gsl_vector_add(SRX,SX);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SX, CH(H02,gsl_matrix_get(C,j,0)));

        gsl_vector_sub(SRX,SX);

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SRX, exp(MulVV(SX,gammai))*gsl_matrix_get(FUNE,1,j));


        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+p2+1,gsl_vector_get(SRX,q));
        }

        if((int)gsl_matrix_get(C,j,1) ==2)
        {
            
            gsl_vector_set_zero(SX);
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                    gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q));
                    gsl_vector_add(SX,RX);
             
                }
            }
            gsl_vector_scale(SX, 1/temp);

            for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));          
            gsl_vector_sub(X,SX); 
             
            gsl_vector_add(X, SRX);

            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+p2+1,gsl_vector_get(X,q));

        }




        /* calculate score for vee */

        double su, uu, sru; 



        /*  vee2 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));

        sru=0;
        for(r=0;r<b;r++)
        {
            if(gsl_matrix_get(H02,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                su=0;
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H02,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                        su+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNUE,0,q);               
                    }
                }

                sru=sru+su*gsl_matrix_get(H02,1,r)/(temp*temp);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        sru=sru*gsl_matrix_get(FUNE,1,j)*exp(MulVV(SX,gammai));

        su=CH(H02,gsl_matrix_get(C,j,0))*gsl_matrix_get(FUNUE,0,j)*exp(MulVV(SX,gammai));

        sru=sru-su;

        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            gsl_vector_set(S,p1+2*p2+1,sru);
        }

        if((int)gsl_matrix_get(C,j,1) ==2)
        {           
            su=0;
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                    su+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNUE,0,q);
                }
            }
            su=su/temp;

            su=gsl_vector_get(FUNU,j)-su;
            su=su+sru;

            gsl_vector_set(S,p1+2*p2+1,su);

        }



        /*  Sigma matrix  */


        for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNBS,t,j));    
        for(i=1;i<p1a;i++)
        {
            for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j));
        }

        for(t=0;t<p1a;t++)
        {
            for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
        }

        for(t=0;t<p1a;t++)
        {
            for(r=0;r<p1a;r++)   gsl_matrix_set(VI,t,r,gsl_matrix_get(bs,t,r));
        }

        gsl_matrix_set(VI,p1a,p1a,gsl_vector_get(FUNUS,j));
        for(t=0;t<p1a;t++)   gsl_matrix_set(VI,t,p1a,gsl_matrix_get(FUNBU,t,j));
        for(t=0;t<p1a;t++)   gsl_matrix_set(VI,p1a,t,gsl_matrix_get(FUNBU,t,j));        


        MulMM(VC,VI,HE);
        MulMM(HE,VC,VI);

        if(p1a==1)
        {
            gsl_vector_set(S,p1+2*p2+2,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+2*p2+3,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+2*p2+4,gsl_matrix_get(VI,1,0)-gsl_matrix_get(VC,1,0));
        }

        if(p1a==2)
        {
            gsl_vector_set(S,p1+2*p2+2,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+2*p2+3,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+2*p2+4,(gsl_matrix_get(VI,2,2)-gsl_matrix_get(VC,2,2))/2);
            gsl_vector_set(S,p1+2*p2+5,gsl_matrix_get(VI,0,1)-gsl_matrix_get(VC,0,1));
            gsl_vector_set(S,p1+2*p2+6,gsl_matrix_get(VI,1,2)-gsl_matrix_get(VC,1,2));
            gsl_vector_set(S,p1+2*p2+7,gsl_matrix_get(VI,0,2)-gsl_matrix_get(VC,0,2));
        }

        if(p1a==3)
        {
            gsl_vector_set(S,p1+2*p2+2,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+2*p2+3,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+2*p2+4,(gsl_matrix_get(VI,2,2)-gsl_matrix_get(VC,2,2))/2);
            gsl_vector_set(S,p1+2*p2+5,(gsl_matrix_get(VI,3,3)-gsl_matrix_get(VC,3,3))/2);
            gsl_vector_set(S,p1+2*p2+6,gsl_matrix_get(VI,0,1)-gsl_matrix_get(VC,0,1));
            gsl_vector_set(S,p1+2*p2+7,gsl_matrix_get(VI,1,2)-gsl_matrix_get(VC,1,2));
            gsl_vector_set(S,p1+2*p2+8,gsl_matrix_get(VI,2,3)-gsl_matrix_get(VC,2,3));
            gsl_vector_set(S,p1+2*p2+9,gsl_matrix_get(VI,0,2)-gsl_matrix_get(VC,0,2));
            gsl_vector_set(S,p1+2*p2+10,gsl_matrix_get(VI,1,3)-gsl_matrix_get(VC,1,3));
            gsl_vector_set(S,p1+2*p2+11,gsl_matrix_get(VI,0,3)-gsl_matrix_get(VC,0,3));
        }



        MulV(S,SS);
        gsl_matrix_add(Cov,SS);
        gsl_vector_add(TS,S);

    }





    gsl_vector_free(FUNU);
    gsl_vector_free(FUNUS);
    gsl_matrix_free(FUNB);
    gsl_matrix_free(FUNBS);
    gsl_matrix_free(FUNBU);
         

    gsl_matrix_free(FUNE); 
    gsl_matrix_free(FUNUSE); 
    gsl_matrix_free(FUNUE); 

    gsl_vector_free(Z);                        
    gsl_vector_free(SZ);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(SZZ);
    gsl_vector_free(X);
    gsl_vector_free(RX);                       
    gsl_vector_free(SX);           
    gsl_vector_free(SRX);
    gsl_vector_free(gammai);

    gsl_vector_free(bi);
    gsl_vector_free(xtilde);
    gsl_vector_free(xtilde1);
    gsl_matrix_free(bs);


    gsl_matrix_free(VC);
    gsl_matrix_free(VI);
    gsl_matrix_free(HE);
    gsl_vector_free(S);                        
    gsl_matrix_free(SS);

    gsl_vector_free(TS);


    return 0;
}




int GetE(
          gsl_vector *FUNU,
          gsl_vector *FUNUS,
          gsl_matrix *FUNB,
          gsl_matrix *FUNBS,
          gsl_matrix *FUNBU,
          gsl_matrix *FUNE,
          gsl_matrix *FUNUSE,
          gsl_matrix *FUNUE,
          gsl_matrix *FUNW,
          gsl_matrix *FUNWB,
          gsl_matrix *FUNWBS,
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl
          )

{
    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;    
    
    int n1 = Y->size1;
    int k = M1->size;

    int i,j,q,t,m,p,status;
    double mu,dem,temp,delta;
    double cuh01,cuh02,haz01,haz02,xgamma1,xgamma2;


    gsl_vector *Z = gsl_vector_alloc(p1),
               *X = gsl_vector_alloc(p2),
               *xtilde = gsl_vector_alloc(p1a),
               *gammai = gsl_vector_alloc(p2);


    gsl_vector_set_zero(FUNU);
    gsl_vector_set_zero(FUNUS);
    gsl_matrix_set_zero(FUNB);
    gsl_matrix_set_zero(FUNBS);
    gsl_matrix_set_zero(FUNBU);
    gsl_matrix_set_zero(FUNE);
    gsl_matrix_set_zero(FUNUSE);
    gsl_matrix_set_zero(FUNUE);
    gsl_matrix_set_zero(FUNW);
    gsl_matrix_set_zero(FUNWB);
    gsl_matrix_set_zero(FUNWBS);



    int point=20;
    int db0,db1,db2,du;


    gsl_vector *xi = gsl_vector_alloc(point);
    gsl_vector *wi = gsl_vector_alloc(point);



    for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

    for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));



    gsl_matrix *VC=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *S=gsl_vector_alloc(p1a+1);
    gsl_matrix *V=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_matrix *VV=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_alloc(p1a+1);


    gsl_matrix_memcpy(VC,sig);


    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }



    gsl_linalg_SV_decomp (VC,V,S,W);
    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++)  gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));
    MulMM(VC,V,VV);


    gsl_matrix_scale(VV,sqrt(2));


    gsl_vector *ri=gsl_vector_alloc(p1a+1);
    gsl_vector *ti=gsl_vector_alloc(p1a+1);


    m=0;    
    for(j=0;j<k;j++)
    {
        dem=0;

        q=(int)gsl_vector_get(M1,j);

        cuh01=CH(H01,gsl_matrix_get(C,j,0));
        cuh02=CH(H02,gsl_matrix_get(C,j,0));

        haz01=HAZ(H01,gsl_matrix_get(C,j,0));
        haz02=HAZ(H02,gsl_matrix_get(C,j,0));


        for(i=0;i<p2;i++) 
        {
            gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
            gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i)); 
        }
        xgamma1=MulVV(X,gammai);

        for(i=0;i<p2;i++) gsl_vector_set(gammai,i,gsl_matrix_get(gamma,1,i));
        xgamma2=MulVV(X,gammai);


        if(p1a==1)
        {

            for(db0=0;db0<point;db0++)
            {
                for(du=0;du<point;du++)
                {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,du));


                    MulM(VV,ri,ti);


                    temp=exp(10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0))));
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));


                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }

                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }


                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+1+p1a));
                        mu=MulVV(Z,beta);
                        
                        delta=gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_matrix_get(Y,m+i,1))/sigma;
                        gsl_matrix_set(FUNW,i,j,gsl_matrix_get(FUNW,i,j)+temp*(kappa+1)/(kappa+delta));


                        for(p=0;p<p1a;p++)  
                        {
                            gsl_matrix_set(FUNWB,i*p1a+p,j,gsl_matrix_get(FUNWB,i*p1a+p,j)+temp*gsl_vector_get(ti,p)*(kappa+1)/(kappa+delta));
                            gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p,j,gsl_matrix_get(FUNWBS,i*p1a*(p1a+1)/2+p,j)+temp*gsl_pow_2(gsl_vector_get(ti,p))
                                           *(kappa+1)/(kappa+delta));
                        }
  
                        for(p=1;p<p1a;p++)
                        {
                            for(t=0;t<p1a-p;t++)   gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j,gsl_matrix_get(FUNWBS,
                                                   i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j)+temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+p)*(kappa+1)/(kappa+delta));
                        }
                    }
                }
            }
        }




        if(p1a==2)
        {

            for(db0=0;db0<point;db0++)
            {
                for(db1=0;db1<point;db1++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=exp(10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0)
                                  -gsl_vector_get(ti,1)*gsl_vector_get(xtilde,1))));
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));
                    
                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));



                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }

                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+1+p1a));
                        mu=MulVV(Z,beta);
                        
                        delta=gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_matrix_get(Y,m+i,1)
                                        -gsl_vector_get(ti,1)*gsl_matrix_get(Y,m+i,2))/sigma;
                        gsl_matrix_set(FUNW,i,j,gsl_matrix_get(FUNW,i,j)+temp*(kappa+1)/(kappa+delta));


                        for(p=0;p<p1a;p++)  
                        {
                            gsl_matrix_set(FUNWB,i*p1a+p,j,gsl_matrix_get(FUNWB,i*p1a+p,j)+temp*gsl_vector_get(ti,p)*(kappa+1)/(kappa+delta));
                            gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p,j,gsl_matrix_get(FUNWBS,i*p1a*(p1a+1)/2+p,j)+temp*gsl_pow_2(gsl_vector_get(ti,p))
                                           *(kappa+1)/(kappa+delta));
                        }
  
                        for(p=1;p<p1a;p++)
                        {
                            for(t=0;t<p1a-p;t++)   gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j,gsl_matrix_get(FUNWBS,
                                                   i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j)+temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+p)*(kappa+1)/(kappa+delta));
                        }

                    }
                }
            }
        }
        }


        if(p1a==3)
        {

        for(db0=0;db0<point;db0++)
        {
            for(db1=0;db1<point;db1++)
            {
                for(db2=0;db2<point;db2++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,db2));
                    gsl_vector_set(ri,3,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=exp(10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0)
                                  -gsl_vector_get(ti,1)*gsl_vector_get(xtilde,1)-gsl_vector_get(ti,2)*gsl_vector_get(xtilde,2))));
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,db2)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));
                    
                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));


                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }

                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }


                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+1+p1a));
                        mu=MulVV(Z,beta);
                        
                        delta=gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_matrix_get(Y,m+i,1)
                                        -gsl_vector_get(ti,1)*gsl_matrix_get(Y,m+i,2)-gsl_vector_get(ti,2)*gsl_matrix_get(Y,m+i,3))/sigma;
                        gsl_matrix_set(FUNW,i,j,gsl_matrix_get(FUNW,i,j)+temp*(kappa+1)/(kappa+delta));


                        for(p=0;p<p1a;p++)  
                        {
                            gsl_matrix_set(FUNWB,i*p1a+p,j,gsl_matrix_get(FUNWB,i*p1a+p,j)+temp*gsl_vector_get(ti,p)*(kappa+1)/(kappa+delta));
                            gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p,j,gsl_matrix_get(FUNWBS,i*p1a*(p1a+1)/2+p,j)+temp*gsl_pow_2(gsl_vector_get(ti,p))
                                           *(kappa+1)/(kappa+delta));
                        }
  
                        for(p=1;p<p1a;p++)
                        {
                            for(t=0;t<p1a-p;t++)   gsl_matrix_set(FUNWBS,i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j,gsl_matrix_get(FUNWBS,
                                                   i*p1a*(p1a+1)/2+p1a+t+(p-1)*(p1a-1),j)+temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+p)*(kappa+1)/(kappa+delta));
                        }
                    }
                }
            }
        }
        }
        }




        if(dem==0) return 100;


        gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)/dem);
        gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)/dem);

        for(i=0;i<p1a;i++)
        {
            gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)/dem);
            gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)/dem);
        }

        for(i=0;i<p1a*(p1a+1)/2;i++)
        {
            gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)/dem);
        }


        for(i=0;i<maxl;i++)
        {
            gsl_matrix_set(FUNW,i,j,gsl_matrix_get(FUNW,i,j)/dem);
        }

        for(i=0;i<maxl*p1a;i++)
        {
            gsl_matrix_set(FUNWB,i,j,gsl_matrix_get(FUNWB,i,j)/dem);
        }

        for(i=0;i<maxl*p1a*(p1a+1)/2;i++)
        {
            gsl_matrix_set(FUNWBS,i,j,gsl_matrix_get(FUNWBS,i,j)/dem);
        }


        gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)/dem);
        gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)/dem);

        gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)/dem);
        gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)/dem);

        m+=q;


    }
 

    gsl_vector_free(Z);
    gsl_vector_free(X);
    gsl_vector_free(xtilde);
    gsl_vector_free(gammai);

    gsl_vector_free(xi);
    gsl_vector_free(ti);
    gsl_vector_free(ri);
    gsl_vector_free(wi);
    gsl_matrix_free(VC);  
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(VV);
    gsl_vector_free(W);  

    return 0;
}



double Getloglik(
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl
          )

{
    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;    
    
    int n1 = Y->size1;
    int k = M1->size;

    int i,j,q,t,m,status;
    double mu,temp1,temp;
    double cuh01,cuh02,haz01,haz02,xgamma1,xgamma2;


    gsl_vector *Z = gsl_vector_alloc(p1),
               *X = gsl_vector_alloc(p2),
               *xtilde = gsl_vector_alloc(p1a),
               *gammai = gsl_vector_alloc(p2);



    int point=20;
    int db0,db1,db2,du;


    gsl_vector *xi = gsl_vector_alloc(point);
    gsl_vector *wi = gsl_vector_alloc(point);


    for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

    for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));



    gsl_matrix *VC=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *S=gsl_vector_alloc(p1a+1);
    gsl_matrix *V=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_matrix *VV=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_alloc(p1a+1);


    gsl_matrix_memcpy(VC,sig);


    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }

    gsl_linalg_SV_decomp (VC,V,S,W);
    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++)  gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));
    MulMM(VC,V,VV);


    gsl_matrix_scale(VV,sqrt(2));


    gsl_vector *ri=gsl_vector_alloc(p1a+1);
    gsl_vector *ti=gsl_vector_alloc(p1a+1);


    double loglik=0;

    m=0;    
    for(j=0;j<k;j++)
    {
        q=(int)gsl_vector_get(M1,j);

        cuh01=CH(H01,gsl_matrix_get(C,j,0));
        cuh02=CH(H02,gsl_matrix_get(C,j,0));

        haz01=HAZ(H01,gsl_matrix_get(C,j,0));
        haz02=HAZ(H02,gsl_matrix_get(C,j,0));


        for(i=0;i<p2;i++) 
        {
            gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
            gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i)); 
        }
        xgamma1=MulVV(X,gammai);

        for(i=0;i<p2;i++) gsl_vector_set(gammai,i,gsl_matrix_get(gamma,1,i));
        xgamma2=MulVV(X,gammai);

        temp1=0;

        if(p1a==1)
        {

            for(db0=0;db0<point;db0++)
            {
                for(du=0;du<point;du++)
                {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=gsl_sf_gamma((kappa+1)/2)/(sqrt(M_PI*kappa*sigma)*gsl_sf_gamma(kappa/2))
                              *exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0))));
                    }


                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    temp1+=temp;

                }
            }
        }


        if(p1a==2)
        {

            for(db0=0;db0<point;db0++)
            {
                for(db1=0;db1<point;db1++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);

                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=gsl_sf_gamma((kappa+1)/2)/(sqrt(M_PI*kappa*sigma)*gsl_sf_gamma(kappa/2))
                              *exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0)
                                  -gsl_vector_get(ti,1)*gsl_vector_get(xtilde,1))));
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,du);

                    temp1+=temp;
                    }
                }
            }
        }

        if(p1a==3)
        {

        for(db0=0;db0<point;db0++)
        {
            for(db1=0;db1<point;db1++)
            {
                for(db2=0;db2<point;db2++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,db2));
                    gsl_vector_set(ri,3,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        
                        temp*=gsl_sf_gamma((kappa+1)/2)/(sqrt(M_PI*kappa*sigma)*gsl_sf_gamma(kappa/2))
                              *exp(0-(kappa+1)/2*log(1+1/(sigma*kappa)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0)
                                  -gsl_vector_get(ti,1)*gsl_vector_get(xtilde,1)-gsl_vector_get(ti,2)*gsl_vector_get(xtilde,2))));
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,db2)*gsl_vector_get(wi,du);

                    temp1+=temp;

                    }
                }
            }
        }

        }


        loglik+=log(temp1);

        m+=q;

    }
 

    gsl_vector_free(Z);
    gsl_vector_free(X);
    gsl_vector_free(xtilde);
    gsl_vector_free(gammai);

    gsl_vector_free(xi);
    gsl_vector_free(ti);
    gsl_vector_free(ri);
    gsl_vector_free(wi);
    gsl_matrix_free(VC);  
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(VV);
    gsl_vector_free(W);  

    return loglik;
}









double DiffV(const gsl_vector *veca, const gsl_vector *vecb)
{
    int k=veca->size;
    int i;
    double diff=0;

    for(i=0;i<k;i++)
    {
        if(Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i))>diff) diff=Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i));
    }

    return (diff);   
}


double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
{
    int nrow=matrixa->size1, ncol=matrixa->size2;
    int i, j;
    double diff=0;

    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol;j++)
        {
            if(Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j))>diff) 
               diff=Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j));
        }
    }

    return (diff);
}



double Abs(const double a, const double b)
{
     
    if (a>=b) return a-b;
    else return b-a;
}



int inv_matrix(gsl_matrix *x_square)
{
   int i,j;
   int k = x_square->size1;

   int status;

   gsl_vector *temp_vector=gsl_vector_alloc(k),
              *solution=gsl_vector_alloc(k);
   gsl_matrix *out = gsl_matrix_alloc(k,k);

   for(i=0;i<k;i++)
   {
       for(j=0;j<k;j++) gsl_matrix_set(out,i,j,gsl_matrix_get(x_square,i,j));
   }

   status=gsl_linalg_cholesky_decompn(out);
   if(status==100) return status;

   for (i = 0; i < k; i++)
   {
       gsl_vector_set_all(temp_vector,0);
       gsl_vector_set(temp_vector,i,1);

       status=gsl_linalg_cholesky_solven(out, temp_vector, solution);
       if(status==100) return status;

       gsl_matrix_set_col(x_square,i,solution);
   }
   
   gsl_vector_free(temp_vector);
   gsl_vector_free(solution);
   gsl_matrix_free(out);

   return 0;

}



void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta)
{
    int p = XX->size1;
    int q = XX->size2;

    int i,j;
    double temp;

    for(i=0;i<p;i++)
    {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
    }

}
         

void MulV(const gsl_vector *Z,gsl_matrix *ZZ)
{
    int p = Z->size;
    int i,j;

    for(i=0;i<p;i++)
    {
        for(j=0;j<p;j++) gsl_matrix_set(ZZ,i,j,gsl_vector_get(Z,i)*gsl_vector_get(Z,j));
    }
}
   

double MulVV(const gsl_vector *Z,const gsl_vector *beta)
{
    int p=Z->size;
    int i;
    double temp=0;

    for(i=0;i<p;i++)  temp+=gsl_vector_get(Z,i)*gsl_vector_get(beta,i);

    return (temp);
}


void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB)
{
    int p=A->size1;
    int q=A->size2;
    int k=B->size2;

    int i,j,t;
    double temp;

    for(i=0;i<p;i++)  
    {
        for(j=0;j<k;j++)
        {
            temp=0;
            for(t=0;t<q;t++)  temp+=gsl_matrix_get(A,i,t)*gsl_matrix_get(B,t,j);
            gsl_matrix_set(AB,i,j,temp);
        }
    }

}


double CH(const gsl_matrix *H, double t)
{
    int a=H->size2;
    int i;

    double ch;

    if(t<gsl_matrix_get(H,0,0)) ch=0;
    else
    {
        ch=0;
        i=0;
        do
        {
            ch+=gsl_matrix_get(H,2,i);
            i+=1;
        }while(i<=a-1 && t>=gsl_matrix_get(H,0,i));
    }

    return (ch);
}


double HAZ(const gsl_matrix *H, double t)
{
    int a=H->size2;
    int i;
    double temp=0;

    for(i=0;i<a;i++)
    {
        if(t==gsl_matrix_get(H,0,i)) temp=gsl_matrix_get(H,2,i);
    }

    return (temp);
}       




double Min(const double t1, const double t2)
{
    if(t1<t2) return t1;
    else return t2;
}


int Diff(
         const gsl_vector *prebeta,
         const gsl_vector *beta,
         const gsl_matrix *pregamma,
         const gsl_matrix *gamma,
         const gsl_vector *prevee,
         const gsl_vector *vee,
         const gsl_matrix *preH01,
         const gsl_matrix *H01,
         const gsl_matrix *preH02,
         const gsl_matrix *H02,
         const double presigma,
         const double sigma,
         const gsl_matrix *presig,
         const gsl_matrix *sig
         )
{    

     double epsilon=0.0001;  
        
     if(DiffV(prebeta,beta)>epsilon || DiffM(pregamma,gamma)>epsilon || DiffV(prevee,vee)>epsilon || DiffM1(preH01,H01)==1 || 
           DiffM1(preH02,H02)==1 || Abs(presigma,sigma)>epsilon  || DiffM(presig,sig)>epsilon)
 
     return 1;

     else return 0;

}


int Sbeta(gsl_vector *beta, double *sigma, const gsl_matrix *Y, const int p1a)
{
    *sigma=0;

    int i,j;

    int p1=Y->size2-1-p1a;
    int n=Y->size1;

    gsl_matrix *X = gsl_matrix_alloc(n,p1);
    gsl_matrix *TX= gsl_matrix_alloc(p1,n);
    gsl_matrix *XX= gsl_matrix_alloc(p1,p1);
    gsl_vector *XY= gsl_vector_alloc(p1);

    gsl_vector *sy = gsl_vector_alloc(n);
    gsl_vector *py = gsl_vector_alloc(n);

    for(i=0;i<n;i++)
    {
        gsl_vector_set(sy,i,gsl_matrix_get(Y,i,0));

        for(j=0;j<p1;j++)
        {
            gsl_matrix_set(X,i,j,gsl_matrix_get(Y,i,p1a+1+j));
        }
    }        
    

    TransM(X,TX);
    MulMM(TX,X,XX);

    int status=inv_matrix(XX);
    if(status==100) return status;

    MulM(TX,sy,XY);
    MulM(XX,XY,beta);
    MulM(X,beta,py);
    gsl_vector_sub(sy,py);

    *sigma=MulVV(sy,sy)/(double)(n-p1);
     

    gsl_matrix_free(X);
    gsl_matrix_free(TX);
    gsl_matrix_free(XX);
    gsl_vector_free(XY);
    gsl_vector_free(sy);
    gsl_vector_free(py);

    return 0;

}







int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
{
    int nrow=matrixa->size1, ncol=matrixa->size2;
    int i, j;
    int diff=0;

    i=nrow-1;

        for(j=0;j<ncol;j++)
        {
            if(gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)>0.001 || gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)<-0.001)

                diff=1;

        }

    return diff;
}


int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1,const gsl_matrix *preH2,const gsl_matrix *H2)
{
     
     if(DiffM1(preH1,H1)==1 || DiffM1(preH2,H2)==1 ) return 1;
     else return 0;
}


void TransM(const gsl_matrix *A, gsl_matrix *B)
{
    int rowa = A->size1;
    int cola = A->size2;

    int i, j;

    for(i=0;i<rowa;i++)
    {
        for(j=0;j<cola;j++)
        {
            gsl_matrix_set(B,j,i,gsl_matrix_get(A,i,j));
        }
    }

}   


void STAT(gsl_matrix *store,int i,double *mean,double *sd)
{
    int n=store->size1;
    int j;

    *mean=0, *sd=0;
    
    for(j=0;j<n;j++)   *mean+=gsl_matrix_get(store,j,i);
    *mean=*mean/(double)n;

    for(j=0;j<n;j++)   *sd+=gsl_pow_2(gsl_matrix_get(store,j,i)-*mean);
    *sd = sqrt(*sd/(double)(n-1));

}
        
     

