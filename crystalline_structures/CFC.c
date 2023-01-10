#include <stdio.h>
#define MAX 100000
int main(){
	int i,j,k;
	float rx[MAX],ry[MAX],rz[MAX];
	FILE *in;
	
	int   cont = 1;
	int   NC   = 10; //número de células unitárias
	float a  = 1.0;
	
	in = fopen("cfc.xyz","w");
	fprintf(in,"%d\n\n",4*NC*NC*NC);
	
		for (i=1; i<=NC;i++){
			for(j=1; j<=NC;j++){
				for(k=1;k<=NC;k++){
					rx[cont]=a*(1.0*(i-1));
					ry[cont]=a*(1.0*(j-1));
					rz[cont]=a*(1.0*(k-1));
					fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
					cont++;
				}
			}
		}
		
		for(i=1; i<=NC;i++){
			for(j=1; j<=NC;j++){
				for(k=1;k<=NC;k++){
					rx[cont]=a*0.5*(1.0*(2*i-1));
					ry[cont]=a*0.5*(1.0*(2*j-1));
					rz[cont]=a*(1.0*(k-1));
					fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
					cont++;
				}
			}
		}

		for(i=1; i<=NC;i++){
			for(j=1; j<=NC;j++){
				for(k=1;k<=NC;k++){
					rx[cont]=a*0.5*(1.0*(2*i-1));
					ry[cont]=a*(1.0*(j-1));
					rz[cont]=a*0.5*(1.0*(2*k-1));
					fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
					cont++;
				}
			}
		}

		for (i=1; i<=NC;i++){
			for(j=1; j<=NC;j++){
				for(k=1;k<=NC;k++){
					rx[cont]=a*(1.0*(i-1));
					ry[cont]=a*0.5*(1.0*(2*j-1));
					rz[cont]=a*0.5*(1.0*(2*k-1));
					fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
					cont++;
				}
			}
		}
		
	fclose(in);
}

