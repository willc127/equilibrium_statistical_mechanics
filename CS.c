#include <stdio.h>
#define MAX 100000
int main(){
	int i,j,k;
	float rx[MAX],ry[MAX],rz[MAX];
	FILE *in;
	
	int cont=1;
	int NC = 10;
	float a =1.0;
	
	//create a file with the coordinates of the molecules in format .xyz
	in=fopen("conf.xyz","w");
	fprintf(in,"%d\n\n",NC*NC*NC);

	//Structure: Simple cubic
	for (i=1; i<=NC;i++){
		for(j=1; j<=NC;j++){
			for(k=1;k<=NC;k++){
				rx[cont]=a*(1.0*(i-1)); //position in x direction
				ry[cont]=a*(1.0*(j-1)); //position in y direction
				rz[cont]=a*(1.0*(k-1)); //position in z direction
				fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
				cont++;
			}
		}
	}
	fclose(in);
}

