#include <stdio.h>
#include <math.h>
#define max 10000

void accel(float m,float kb,int n,float eps24,float eps,float rcutsq,
		float rx[max],float ry[max],float rz[max],float ax[max],float ay[max],float az[max],
		float sigsq, float boxlx, float boxly, float boxlz, float *V, float *W);
		
float potential(float rijsq,float sigsq, float *wij, float eps24, float kb);

int main(){
	int i,n,k;
	float rx[max],ry[max],rz[max],vx[max],vy[max],vz[max],ax[max],ay[max],az[max];
	char dummy; //para ler a linha vazia do arquivo cfc
	FILE *in;
	
	int atom = 12; //part�cula a ser utilizada
	n = 864;  // n�mero de part�culas
	in=fopen("cfc.xyz","r"); //usar'r' somente para leitura
	fscanf(in,"%d",&n);
	printf("%d\n",n);
	fscanf(in,"%c",&dummy);
	
	for(i=0;i<n;i++){		
	fscanf(in,"%d %f %f %f %f %f %f",&atom,&rx[i],&ry[i],&rz[i],&vx[i],&vy[i],&vz[i]); //o caractere & serve para salvar o valor da vari�vel a cada rodada
}
fclose(in);

	float m   = 16.0; //massa molar
	float NA = 6.0221409E23;
	float dens = 100.0;// Densidade em kg/m�
	float a = pow(4.0*m/(1000.0*dens*NA),1.0/3.0)*1e10; //aresta em angstrom
	float kb  = 1.38064852E-23;
	float r   = NA*kb;
	float sigma = 3.73; //angstrom
	float eps = 148.0; // em K
	float eps24 = 24*eps*r/(m*1e-3);
	float rcutsq = 2.5*sigma*2.5*sigma;
	int NC = cbrt(n*0.25); //n�mero de c�lulas
	printf("nc = %d\n",NC);
	float sigsq=sigma*sigma;
	float boxlx = a*NC; //define a aresta total do cubo
	float boxly = a*NC;
	float boxlz = a*NC;
	float V ;
	float W;
	float wij;

	accel(m,kb,n,eps24,eps,rcutsq,rx,ry,rz,ax,ay,az,sigsq,boxlx,boxly,boxlz,&V,&W);

	FILE *out;
	out = fopen("movs.xyz","w");
	
	float t = 0;
	float dt = 1;
	float dt2 = dt*0.5;
	float tf = 6000;
	float T;
	float kin;
	float pot;
	int ct = 0;
	FILE *temp;
	temp = fopen("temper.dat","w");

	while(t<tf){	
		for (i=0;i<n;i++){
			vx[i]=vx[i]+ax[i]*dt2;
			vy[i]=vy[i]+ay[i]*dt2;
			vz[i]=vz[i]+az[i]*dt2;
			
			rx[i]=rx[i]+vx[i]*dt;
			ry[i]=ry[i]+vy[i]*dt;
			rz[i]=rz[i]+vz[i]*dt;
		}
			
		accel(m,kb,n,eps24,eps,rcutsq,rx,ry,rz,ax,ay,az,sigsq,boxlx,boxly,boxlz,&V,&W);
			
		float sumv2 = 0.0;
	
		for (i=0;i<n;i++){
			vx[i]=vx[i]+ax[i]*dt2;
			vy[i]=vy[i]+ay[i]*dt2;
			vz[i]=vz[i]+az[i]*dt2;
			sumv2 = sumv2+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];

			if(ct==2){
							
				fprintf(out,"%d\n\n",n);
					
				for (k=0;k<n;k++){ 
					fprintf(out,"12 %f %f %f  \n",rx[k],ry[k],rz[k]);
				//	printf("12 %f %f %f %e %e %e %e %e %e \n",rx[i],ry[i],rz[i],vx[i],vy[i],vz[i],ax[i],ay[i],az[i]);
					ct = 0;	
				}
			}	
		}
	
	ct = ct+1;
	t  = t+dt;
	float rijsq;
	float sigsq;
	float kin = 0.5*m*sumv2*1e7/(float)(n); //J/molecule
	float H = kin+V;
	float T = 2.0*kin/(3.0*r);
	float Vol = a*NC*a*NC*a*NC*1e-30; //m3
	float Pr = ((float)(n)/NA*r*T/Vol + W/3*eps24*m*1e-3/Vol/NA)*1e-5; //bar
	printf("%f %f %f %f %f %f\n",t,T,Pr,kin,V,H);
	fprintf(temp,"%f %f %f %f %f %f\n",t,T,Pr,kin,V,H);
	}
	
fclose(out);
return 0;	
}
	
void accel(float m,float kb,int n, float eps24,float eps,float rcutsq,
	float rx[max],float ry[max],float rz[max],float ax[max],float ay[max],float az[max],
	float sigsq, float boxlx, float boxly, float boxlz, float *V, float *W){

	float rxi,ryi,rzi,aij,rxij,ryij,rzij,pot,rijsq,wij;
	int i,j;
	
	for(i=0;i<n;i++){
		ax[i]=0.0;
		ay[i]=0.0;
		az[i]=0.0;
	}	

	*V = 0.0;
	*W = 0.0;
	for(i=0;i<n-1;i++){
		rxi = rx[i];
		ryi = ry[i];
		rzi = rz[i];
		for(j=i+1;j<n;j++){
			rxij = rxi - rx[j]; 
			ryij = ryi - ry[j]; 
			rzij = rzi - rz[j];
			rxij = rxij-boxlx*round(rxij/boxlx);
			ryij = ryij-boxly*round(ryij/boxly);
			rzij = rzij-boxlz*round(rzij/boxlz);
			rijsq = rxij*rxij+ryij*ryij+rzij*rzij;
				
				if(rijsq<=rcutsq){ 
					pot = potential(rijsq,sigsq,&wij,eps,kb);
					*V = *V+pot;
					aij = wij/rijsq*1e-10; // A/fs�
					ax[i]= ax[i]+aij*rxij;
					ay[i]= ay[i]+aij*ryij;
					az[i]= az[i]+aij*rzij;
					ax[j]= ax[j]-aij*rxij;
					ay[j]= ay[j]-aij*ryij;
					az[j]= az[j]-aij*rzij;	
					*W = *W + wij;
				}
			}
		}
		*V = *V*eps24*m*1e-3/6.0/(float)(n);	
		for (i=0;i<n;i++){
			ax[i]=ax[i]*eps24;
			ay[i]=ay[i]*eps24;
			az[i]=az[i]*eps24;
		}
}

float potential(float rijsq,float sigsq, float *wij,float eps, float kb){
		float sr2 = sigsq/rijsq;
		float sr6 = sr2*sr2*sr2;
		float sr12 = sr6*sr6;
		float pot = sr12-sr6;

		*wij = 2.0*sr12-sr6;
		
		return pot;
	}

