// Aluno Rodrigo Amaral Coutinho Bartolomeu. RA. 209688
// ******************************************************************
// Calculo da G(r)
// ******************************************************************
// Prof. Dr. Luis Fernando Mercier Franco
// Universidade Estadual de Campinas, UNICAMP. 
// Faculdade de Engenharia Quimica                          11/11/2017
// *******************************************************************

#include <stdio.h>
#include <math.h>
#include <string.h>
#define MAX 10000

int main (){
   	int i,j,k,atom,step;
    	float rx[MAX],ry[MAX],rz[MAX];
    	float vx,vy,vz;
	char dummy[100];
   	int hist[MAX];
	float Gr[MAX];
    	float rxij,ryij,rzij,rijsq,rij,x,y,z;
	int n		= 1372;
	int nstep 	= 4001;
	int nstep_equil = 2500;
	float m       	= 0.016;                  // Massa Molecular do Metano (kg.mol⁻¹)
	float na      	= 6.0221409e23;           // Número de Avogadro (mol⁻¹)'	
	float rho   	= 100;                  // Densidade em kg.m⁻³
	float V   = (n*m)/(rho*na);         // Volume da caixa (m³)
	float boxl  	= pow((float)(V),(float)(1.0/3.0))*1e10; // lado da caixa (A)
	float hbox   	= 0.5*boxl;            // Metade do lado da caixa (A)
	float mols    	= na/((float)(n));
	int bin;
    	int maxbin 	= 1000;
    	float delr 	= hbox/(float)(maxbin);
	float constante	= 4.0*M_PI*rho/3.0;
	float rlower,rupper,nideal;

    	FILE *inp;
    	inp = fopen("movs.xyz","r");
   	FILE *in3;
	in3 = fopen("Gr.dat","w");	
	
 	for (bin=0;bin<maxbin;bin++);{
            		hist[bin] = 0;
        	}

	for (step=0;step<nstep;step++){
		printf("Step %d ...\n",step);
		fscanf(inp,"%d",&n);
		fscanf(inp,"%s",dummy);
		for (i=0;i<n;i++){
			fscanf(inp,"%d %f %f %f",&atom,&x,&y,&z);
			if(step>nstep_equil){
		        	rx[i] = x;
        	       		ry[i] = y;
                		rz[i] = z;
			}
	    
	    }
       
        	if(step>nstep_equil){
            		for (i=0;i<n-1;i++){
                		for (j=i+1;j<n;j++){
                    			rxij = rx[i] - rx[j];
		                	ryij = ry[i] - ry[j];
                			rzij = rz[i] - rz[j];
                    			rxij = rxij -boxl*round(rxij/boxl);
                    			ryij = ryij -boxl*round(ryij/boxl);
                    			rzij = rzij -boxl*round(rzij/boxl);
                    			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
                    			rij = sqrt(rijsq);
                    			bin = (int)(rij/delr);
                    			if (bin<=maxbin){
                        			hist[bin]= hist[bin] + 2;
                   			}
			
               		 	}	        
        		}
			for (bin=0;bin<maxbin;bin++){
				rlower  = (float)(bin)*delr;
				rupper  = rlower+delr;
				nideal  = constante*(rupper*rupper*rupper - rlower*rlower*rlower);
				Gr[bin] = (float)(hist[bin])/(float)(nstep-nstep_equil)/(float)(n)/nideal;
				fprintf(in3,"%f %f \n",(float)(bin)*delr,Gr[bin]);

			}
        	}
   			
	}




}

