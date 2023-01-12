#include <stdio.h>
#include <math.h>
#define max 10000

int main(){
	int i,v;
	float t,tf,T[max],kin[max],V[max],H[max],soma=0,somavar=0,somaT=0,somaH =0,somacv=0,somaP=0,Pr[max],Cv;
	FILE *in;
	FILE *intf;
	
	int n = 500;
	float kb  = 1.38064852E-23;
	float NA = 6.0221409E23;
	float dens = 100;
	float m = 16;
	
	in = fopen("temper.dat","r"); //usar'r' somente para leitura
	intf = fopen("tf.dat","r");
	v = 2000; //valor de t onde acaba a zona de equilibra��o
	fscanf(intf,"%f",&tf);
	//printf("%f\n",tf);

	for(i=1;i<=tf;i++){ //para o c�lculo da energia cin�tica m�dia
		
		fscanf(in,"%f %f %f %f %f %f ",&t,&T[i],&Pr[i],&kin[i],&V[i],&H[i]);	
		if (i>=v){
		
		soma = soma + kin[i];
		somaH = somaH + H[i];
		}
	}

		fclose(in);
	
	float kinmedia = soma/(tf-v+1);
	
	float sigma = 373e-10;
	float Dens_corr = dens*sigma*sigma*sigma;
	float rcut_corr = 2.5;
	float pi = 3.14159265;
	float eps  = 148;
	
	float V_corr1 = ((8*pi*n*Dens_corr)/(9*pow(rcut_corr,9.0)))-((8*pi*n*Dens_corr)/(3*pow(rcut_corr,3.0)));
	float V_corr2 = eps/m*1000*V_corr1;
	
	float P_corr1 = (32/9*pi*Dens_corr*Dens_corr)/(pow(rcut_corr,9.0))-((16/3*pi*Dens_corr*Dens_corr)/pow(rcut_corr,3.0))*NA*1e-5;
	float P_corr2 = eps/m*P_corr1/(sigma*sigma*sigma);
		
	for(i=1;i<=tf;i++){ //calcula a vari�ncia da energia cin�tica
		
		if (i>=v){
		somavar = somavar + (kin[i]-kinmedia)*(kin[i]-kinmedia);
		somaP = somaP + Pr[i];
	
		}
	}
	
	float var = somavar/(tf-v+1);

		for(i=1;i<=tf;i++){ //Calcula o valor de cv
					
			if (i>=v){
						
			Cv = (NA*((3*kb)/(2*(1-(var*(2/3)*(1/(n*kb*kb*T[i]*T[i]))))))); //J/molK
			printf("Cv: %f \n", Cv);	
			somacv = somacv + Cv;
			somaT = somaT + T[i];
			
			}
		}
	
	float cvs = somacv/(tf-v+1);
	float Tm = somaT/(tf-v+1);
	float Prm = somaP/(tf-v+1)+P_corr2;
	float Hm = somaH/(tf-v+1);
	
   	float T_cv_v1 = 4170;
   	float T_cv_v2 = 2180;
   	float T_cv_v3 = 4320;
   	float T_cv_v4 = 1870;

    float termo1_1 = (T_cv_v1/Tm);
    float termo1_2 = (T_cv_v2/Tm);
    float termo1_3 = (T_cv_v3/Tm);
    float termo1_4 = (T_cv_v4/Tm);
    
    float den1 = (exp(termo1_1)-1)*(exp(termo1_1)-1);
    float den2 = (exp(termo1_2)-1)*(exp(termo1_2)-1);
    float den3 = (exp(termo1_3)-1)*(exp(termo1_3)-1);
    float den4 = (exp(termo1_4)-1)*(exp(termo1_4)-1);
    
    float termo2_1 = exp(termo1_1)/den1;
    float termo2_2 = exp(termo1_2)/den2;
    float termo2_3 = exp(termo1_3)/den3;
    float termo2_4 = exp(termo1_4)/den4;
    
    float deg1 = 1;
    float deg2 = 2;
    float deg3 = 3;
    float deg4 = 3;
    
    float cv_v = NA*kb*(deg1*(termo1_1*termo1_1)*termo2_1+deg2*(termo1_2*termo1_2)*termo2_2+deg3*(termo1_3*termo1_3)*termo2_3+deg4*(termo1_4*termo1_4)*termo2_4);

	// C�lculo de Cv
    float Cv3   = Cv + cv_v + 1.5*NA*kb;	
	//printf("%f %f %f %f \n",Tm,Prm,Hm+V_corr2,Cv3);	

return 0;	
}


