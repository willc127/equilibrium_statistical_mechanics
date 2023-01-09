#include <stdio.h>
#include <math.h>
#include <string.h>
#define max 10000

int main(){
	
int nstep=4001,nstep_equil=2500,n,nc,i,j,atom,step;
int maxbin=1000;
char dummy;  

float acalc,a,boxl, dens,m,na,x,y,z;
float rx[max], ry[max], rz[max],gr[max];
float rlower,rupper,nIdeal,rij,rxij,ryij,rzij,rijsq;
int hist[max];

FILE *in;
in=fopen("movs.xyz","r");

FILE *out;

na      = 6.0221409e23;
m		=	0.016;
dens 	= 	100.0;
nc		= 	pow(n/4,1.0/3.0);              
acalc	=	(n*m/na)/(dens*nc*nc*nc);  
a		=	(pow(acalc,1.0/3.0))*(1e10);
boxl	= 	nc*a;
     
float delR=0.5*boxl/(float)(maxbin);  
float vol=pow(boxl,3.0/1.0);
float densn=n/vol;
float constante=4*M_PI*densn/3.0;
int bin;
	
	for (bin=0;bin<maxbin;bin++)
		hist[bin]=0;

	for (step=0;step<nstep;step++){
		printf("Step %d ...\n",step);
		fscanf(in,"%d",&n);
		fscanf(in,"%c",&dummy); // %s ou %c?
		for (i=0;i<n;i++)
			fscanf(in, "%d %f %f %f", &atom, &rx[i], &ry[i],&rz[i]);
	
		if(step>=nstep_equil){
			for(i=0;i<n-1;i++){
			
				for(j=i+1;j<n;j++){
				
					rxij 	=rx[i]-rx[j];
					ryij 	=ry[i]-ry[j];
					rzij 	=rz[i]-rz[j];
           	   		rxij 	=rxij-boxl*round(rxij/boxl);
            		ryij 	=ryij-boxl*round(ryij/boxl);
           			rzij 	=rzij-boxl*round(rzij/boxl);
            		rijsq	=rxij*rxij+ryij*ryij+rzij*rzij;
            		rij		=sqrt(rijsq);
            		bin		=(int)(rij/delR);
            		if(bin<=maxbin)
            			hist[bin]=hist[bin]+2;
				}	
			}
		}
	}
	fclose(in);	
	
	out=fopen("Gr.dat","w"); 	
	for (bin=0;bin<=maxbin;bin++){
		rlower	=(float)(bin)*delR;
		rupper	=rlower+delR;
		nIdeal	=constante*((rupper*rupper*rupper)-(rlower*rlower*rlower));
		gr[bin]	=(float)(hist[bin])/(float)(nstep-nstep_equil)/(float)(n)/nIdeal;
		printf("%e %e %e\n",bin,delR,bin*delR);
		//fprintf(out, "%f %e\n",(float)(bin)*delR,gr[bin]);
	}
	fclose(out);
	return 0;	
}
