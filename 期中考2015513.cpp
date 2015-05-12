#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int FFTV3(double *x_r, double *x_i, double *y_r, double *y_i,int N);
int FFTV2(double *x_r, double *x_i, double *y_r, double *y_i,int N);
int FFTV5(double *x_r, double *x_i, double *y_r, double *y_i,int N);

int main()
{
	
	int k,n,N;
	double *x_r, *x_i, *y_r, *y_i, w_r, w_i;
	clock_t t1, t2;
		
	printf("please input N=");
	scanf("%d",&N);
	printf("N=%d\n",N);
		
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));	
	y_i = (double *) malloc(N*sizeof(double));
	for(n=0;n<N;++n)	
	{
		x_r[n] = n;
		x_i[n] = 0;
	}

    t1 = clock();
	if (N % 3 == 0 )
	{
    FFTV3(x_r, x_i, y_r, y_i, N);
    }
	else if (N % 2 ==0)
	{
    FFTV2(x_r, x_i, y_r, y_i, N);
    }
    else 
    {
    FFTV5(x_r, x_i, y_r, y_i, N);
	}
     t2 = clock();
	printf("%f secs\n",1.0*(t2-t1)/CLOCKS_PER_SEC);
	system("pause"); 
	for(n=0;n<N;++n)
	{
		printf("%d : %f + %f i\n",n,y_r[n],y_i[n]);
	}

	return 0;		
}

int FFTV3(double *x_r, double *x_i, double *y_r, double *y_i,int N)
{
	
	if(N==1)
	{
		y_r[0]=x_r[0];
		y_i[0]=x_i[0];
		return 0;
	}
	int k,n;
	double *u_r, *u_i, *v_r, *v_i, w_r, w_i,w2r,w2i,w_a,w_b;
	
	u_r = (double *)malloc(N*sizeof(double));
	u_i = (double *)malloc(N*sizeof(double));
	v_r = (double *)malloc(N*sizeof(double));
	v_i = (double *)malloc(N*sizeof(double));
	for(n=0;n<N/3;n++)
	{
		u_r[n] = x_r[3*n];
		u_i[n] = x_i[3*n];
		u_r[n+N/3] = x_r[3*n+1];
		u_i[n+N/3] = x_i[3*n+1];
		u_r[n+2*N/3] = x_r[3*n+2];
		u_i[n+2*N/3] = x_i[3*n+2];
	}
	if (n % 2 == 0)
    {
       FFTV2(u_r, u_i, v_r, v_i, N/3);
	   FFTV2(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
	   FFTV2(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
       }
       

    else if ((n % 3) == 0) 
    	{
       FFTV3(u_r, u_i, v_r, v_i, N/3);
	   FFTV3(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
	   FFTV3(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
        }
   
    else 
	{
	   FFTV5(u_r, u_i, v_r, v_i, N/3);
	   FFTV5(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
	   FFTV5(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
	}
	  
    for(k=0;k<N/3;++k)
	{   
		w_r=cos(-k*2*M_PI/N);
		w_i=sin(-k*2*M_PI/N);
		w2r=cos(-2*k*2*M_PI/N);
		w2i=sin(-2*k*2*M_PI/N);
		w_a=cos(-2*M_PI/3);
		w_b=sin(-2*M_PI/3);
		//printf("%f %f %f %f %d %d\n",w_r,w_i,w2r,w2i,k,N);
		y_r[k]=v_r[k]+w_r*v_r[k+N/3]-w_i*v_i[k+N/3]+w2r*v_r[k+2*N/3]-w2i*v_i[k+2*N/3];
		
		y_i[k]=v_i[k]+w_r*v_i[k+N/3]+w_i*v_r[k+N/3]+w2r*v_i[k+2*N/3]+w2i*v_r[k+2*N/3];
		y_r[k+N/3]=v_r[k]
		           +(w_r*v_r[k+N/3]-w_i*v_i[k+N/3])*w_a-(w_r*v_i[k+N/3]+w_i*v_r[k+N/3])*w_b
		           +(w2r*v_r[k+2*N/3]-w2i*v_i[k+2*N/3])*cos(-4*M_PI/3)-(w2r*v_i[k+2*N/3]+w2i*v_r[k+2*N/3])*sin(-4*M_PI/3);
		y_i[k+N/3]=v_i[k]
		          +(w_r*v_r[k+N/3]-w_i*v_i[k+N/3])*w_b+(w_r*v_i[k+N/3]+w_i*v_r[k+N/3])*w_a
				  +(w2r*v_r[k+2*N/3]-w2i*v_i[k+2*N/3])*sin(-4*M_PI/3)+(w2r*v_i[k+2*N/3]+w2i*v_r[k+2*N/3])*cos(-4*M_PI/3);
		y_r[k+2*N/3]=v_r[k]
		           +(w_r*v_r[k+N/3]-w_i*v_i[k+N/3])*cos(-4*M_PI/3)-(w_r*v_i[k+N/3]+w_i*v_r[k+N/3])*sin(-4*M_PI/3)
		           +(w2r*v_r[k+2*N/3]-w2i*v_i[k+2*N/3])*w_a-(w2r*v_i[k+2*N/3]+w2i*v_r[k+2*N/3])*w_b;
				   
		y_i[k+2*N/3]=v_i[k]
		            +(w_r*v_r[k+N/3]-w_i*v_i[k+N/3])*sin(-4*M_PI/3)+(w_r*v_i[k+N/3]+w_i*v_r[k+N/3])*cos(-4*M_PI/3)
				    +(w2r*v_r[k+2*N/3]-w2i*v_i[k+2*N/3])*w_b+(w2r*v_i[k+2*N/3]+w2i*v_r[k+2*N/3])*w_a;
	}
	    return 0;
}

int FFTV2(double *x_r, double *x_i, double *y_r, double *y_i,int N)
{
	if(N==1)
	{
		y_r[0]=x_r[0];
		y_i[0]=x_i[0];
		return 0;
	}
	int k,n;
	double *u_r, *u_i, *v_r, *v_i, w_r, w_i;
	
	u_r = (double *)malloc(N*sizeof(double));
	u_i = (double *)malloc(N*sizeof(double));
	v_r = (double *)malloc(N*sizeof(double));
	v_i = (double *)malloc(N*sizeof(double));
	
	for(n=0;n<N/2;++n)
	{
		u_r[n] = x_r[2*n];
		u_i[n] = x_i[2*n];
		u_r[n+N/2] = x_r[2*n+1];
		u_i[n+N/2] = x_i[2*n+1];
	}
	if (n % 2 == 0)
	{
		FFTV2(u_r, u_i, v_r, v_i, N/2);
	    FFTV2(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
	}
	
	else if (n % 3 == 0)
	{
		FFTV3(u_r, u_i, v_r, v_i, N/2);
	    FFTV3(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
	}
	
	else
	{
		FFTV5(u_r, u_i, v_r, v_i, N/2);
	    FFTV5(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
	}
	
	
	for(k=0;k<N/2;++k)
	{
		w_r=cos(-k*2*M_PI/N);
		w_i=sin(-k*2*M_PI/N);
		//printf("%f %f\n",w_r,w_i);
		y_r[k]=v_r[k]+w_r*v_r[k+N/2]-w_i*v_i[k+N/2];
		y_i[k]=v_i[k]+w_r*v_i[k+N/2]+w_i*v_r[k+N/2];
		y_r[k+N/2]=v_r[k]-(w_r*v_r[k+N/2]-w_i*v_i[k+N/2]);
		y_i[k+N/2]=v_i[k]-(w_r*v_i[k+N/2]+w_i*v_r[k+N/2]);
	}
	return 0;
}

int FFTV5(double *x_r, double *x_i, double *y_r, double *y_i,int N)
{

	if(N==1)
	{
		y_r[0]=x_r[0];
		y_i[0]=x_i[0];
		return 0;
	}
	int k,n;
	double *u_r, *u_i, *v_r, *v_i, w_r, w_i,w2r, w2i,w3r,w3i,w4r,w4i,w_a,w_b,w2a,w2b,w3a,w3b,w4a,w4b;
	
	
	
	u_r = (double *)malloc(N*sizeof(double));
	u_i = (double *)malloc(N*sizeof(double));
	v_r = (double *)malloc(N*sizeof(double));
	v_i = (double *)malloc(N*sizeof(double));
    
	for(n=0;n<N/5;n++)
	{
		u_r[n] = x_r[5*n];
		u_i[n] = x_i[5*n];
		u_r[n+N/5] = x_r[5*n+1];
		u_i[n+N/5] = x_i[5*n+1];
		u_r[n+2*N/5] = x_r[5*n+2];
		u_i[n+2*N/5] = x_i[5*n+2];
	    u_r[n+3*N/5] = x_r[5*n+3];
	    u_i[n+3*N/5] = x_i[5*n+3];
		u_r[n+4*N/5] = x_r[5*n+4];
		u_i[n+4*N/5] = x_i[5*n+4];
	    
	}
	
	if (n % 2 == 0)
	{
		FFTV2(u_r, u_i, v_r, v_i, N/5);
        FFTV2(u_r+N/5, u_i+N/5, v_r+N/5, v_i+N/5, N/5);
        FFTV2(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
	    FFTV2(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
	    FFTV2(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}
	
	else if (n % 3 == 0)
	{
		FFTV3(u_r, u_i, v_r, v_i, N/5);
        FFTV3(u_r+N/5, u_i+N/5, v_r+N/5, v_i+N/5, N/5);
        FFTV3(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
    	FFTV3(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
	    FFTV3(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}
	
	else
	{
    	FFTV5(u_r, u_i, v_r, v_i, N/5);
        FFTV5(u_r+N/5, u_i+N/5, v_r+N/5, v_i+N/5, N/5);
        FFTV5(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
	    FFTV5(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
	    FFTV5(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}
	

	for(k=0;k<N/5;++k)
	{
        w_r=cos(-k*2*M_PI/N);
		w_i=sin(-k*2*M_PI/N);
		w2r=cos(-2*k*2*M_PI/N);
		w2i=sin(-2*k*2*M_PI/N);
		w3r=cos(-3*k*2*M_PI/N);
		w3i=sin(-3*k*2*M_PI/N);
		w4r=cos(-4*k*2*M_PI/N);
		w4i=sin(-4*k*2*M_PI/N);
        w_a=cos(-2*M_PI/5);
        w_b=sin(-2*M_PI/5);
        w2a=cos(-2*2*M_PI/5);
        w2b=sin(-2*2*M_PI/5);
        w3a=cos(-3*2*M_PI/5);
        w3b=sin(-3*2*M_PI/5);
        w4a=cos(-4*2*M_PI/5);
        w4b=sin(-4*2*M_PI/5);
        
		y_r[k]=v_r[k]
		       +w_r*v_r[k+N/5]-w_i*v_i[k+N/5]
			   +w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5]
               +w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5]
	           +w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5];
	    y_i[k]=v_i[k]
		       +w_r*v_i[k+N/5]+w_i*v_r[k+N/5]
			   +w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5]
	           +w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5]
	           +w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5];
	    
		y_r[k+N/5]=v_r[k]
		           +w_a*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
                   -w_b*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])	
				   +w2a*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])          
				   -w2b*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])   
				   +w3a*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				   -w3b*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
	               +w4a*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5])
	               -w4b*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5]);
	    
		y_i[k+N/5]=v_i[k]
		           +w_a*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])  
		           +w_b*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
		           +w2a*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])
		           +w2b*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])
	               +w3a*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
                   +w3b*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				   +w4a*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5])	
	               +w4b*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5]);
	   
	    y_r[k+2*N/5]=v_r[k]
		            +w2a*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
                    -w2b*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])	
				    +w4a*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])          
				    -w4b*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])   
				    +w_a*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    -w_b*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
	                +w3a*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5])
	                -w3b*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5]);             
	               
	    y_i[k+2*N/5]=v_i[k]
	                +w2a*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])  
		            +w2b*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
		            +w4a*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])
		            +w4b*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])
	                +w_a*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
                    +w_b*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    +w3a*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5])	
	                +w3b*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5]);
	               
	    y_r[k+3*N/5]=v_r[k]
		            +w3a*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
                    -w3b*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])	
				    +w_a*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])          
				    -w_b*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])   
				    +w4a*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    -w4b*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
	                +w2a*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5])
	                -w2b*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5]);
	               
	    y_i[k+3*N/5]=v_i[k]
	                +w3a*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])  
		            +w3b*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
		            +w_a*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])
		            +w_b*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])
	                +w4a*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
                    +w4b*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    +w2a*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5])	
	                +w2b*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5]);           
	    
		y_r[k+4*N/5]=v_r[k]
		            +w4a*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
                    -w4b*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])	
				    +w3a*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])          
				    -w3b*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])   
				    +w2a*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    -w2b*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
	                +w_a*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5])
	                -w_b*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5]);           
	            
	    y_i[k+4*N/5]=v_i[k]
	                +w4a*(w_r*v_i[k+N/5]+w_i*v_r[k+N/5])  
		            +w4b*(w_r*v_r[k+N/5]-w_i*v_i[k+N/5])
		            +w3a*(w2r*v_i[k+2*N/5]+w2i*v_r[k+2*N/5])
		            +w3b*(w2r*v_r[k+2*N/5]-w2i*v_i[k+2*N/5])
	                +w2a*(w3r*v_i[k+3*N/5]+w3i*v_r[k+3*N/5])
                    +w2b*(w3r*v_r[k+3*N/5]-w3i*v_i[k+3*N/5])
				    +w_a*(w4r*v_i[k+4*N/5]+w4i*v_r[k+4*N/5])	
	                +w_b*(w4r*v_r[k+4*N/5]-w4i*v_i[k+4*N/5]);            
	}
	   return 0;
}

	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	               
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
