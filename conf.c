// *********************************************************************
// Universidade Estadual de Campinas
// Faculdade de Engenharia Quimica
// IQ620 - Termodinamica II
// William Cesar de Oliveira Ribeiro
// RA: 209656
// *********************************************************************
// Atom position generator in a cubic centered face configuration with
// random initial velocities following a Gaussian distribution
// *********************************************************************

#include <stdio.h>
#include <math.h>
#define MAX 100000

float gauss(int seed[1]);
float ranf(int seed[1]);

int main()
{
	// Define the variables
	int i, j, k;
	float rx[MAX], ry[MAX], rz[MAX];

	//opening file to write the positions
	FILE *in;
	in = fopen("cfc.xyz", "w");

	//Define the total number of particles
	int n = 864; // n = (4, 32, 108, 256, 500, 864, 1372, 2048, 2916, 4000...)
	int NC = (int)(cbrt(0.25 * (float)(n))); // total number of cells
	printf("n = %d\n", n);
	printf("nc = %d\n", NC);
	fprintf(in, "%d\n\n", n);

	//Define the chemical and thermodynamical parameters
	float T  = 500.0;		   // temperature in Kelvin
	float kb = 1.38064852E-23; // Boltzmann's constant
	float NA = 6.0221409E23;   // Avogadro's number
	float R  = kb * NA;		   // universal gas constant
	float m  = 16.0;		   // molecular weight of methane (g/mol)
	float dens = 100.0;		   // density in kg/m3
	float csi  = sqrt(3.0*(float)(n)*R*T/(m * 1e-3));
	float a  = pow(4*m*1e-3/(dens*NA),1.0/3.0)*1e10; // cubic edge in angstrom
	
	// Atoms positions
	int cont = 0;
	// ! This loop calculates the positions in the vertices
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < NC; j++)
		{
			for (k = 0; k < NC; k++)
			{
				rx[cont] = a * (float)(i);
				ry[cont] = a * (float)(j);
				rz[cont] = a * (float)(k);
				// fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
				cont++;
			}
		}
	}

	// ! This loop calculates the positions in faces of the XY plane
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < NC; j++)
		{
			for (k = 0; k < NC; k++)
			{
				rx[cont] = a * (0.5 + (float)(i));
				ry[cont] = a * (0.5 + (float)(j));
				rz[cont] = a * (float)(k);
				// fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
				cont++;
			}
		}
	}

	// ! This loop calculates the positions in faces of the XZ plane
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < NC; j++)
		{
			for (k = 0; k < NC; k++)
			{
				rx[cont] = a * (0.5 + (float)(i));
				ry[cont] = a * (float)(j);
				rz[cont] = a * (0.5 + (float)(k));
				// fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
				cont++;
			}
		}
	}

	// ! This loop calculates the positions in faces of the YZ plane
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < NC; j++)
		{
			for (k = 0; k < NC; k++)
			{
				rx[cont] = a * (float)(i);
				ry[cont] = a * (0.5 + (float)(j));
				rz[cont] = a * (0.5 + (float)(k));
				// fprintf(in,"C %f %f %f \n",rx[cont],ry[cont],rz[cont]);
				cont++;
			}
		}
	}
	printf("cont = %d\n", cont);

	// initial velocities calculation
	int seed[1];
	float vx[50000], vy[50000], vz[50000];

	// define a seed for the random number generator
	seed[0] = 408065;

	// calculate the velocities
	float sum = 0.0;
	for (i = 0; i < n; i++)
	{
		vx[i] = gauss(seed);
		vy[i] = gauss(seed);
		vz[i] = gauss(seed);
		sum = sum + vx[i] + vy[i] + vz[i];
	}

	for (i = 0; i < n; i++)
	{
		vx[i] = vx[i] / sum;
		vy[i] = vy[i] / sum;
		vz[i] = vz[i] / sum;
	}

	float sumvsq = 0.0;
	for (i = 0; i < n; i++)
	{
		sumvsq = sumvsq + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	}

	float sumv = sqrt(sumvsq);
	float sumvx = 0.0;
	float sumvy = 0.0;
	float sumvz = 0.0;

	for (i = 0; i < n; i++)
	{
		vx[i] = csi * vx[i] / sumv;
		vy[i] = csi * vy[i] / sumv;
		vz[i] = csi * vz[i] / sumv;
		sumvx = sumvx + vx[i];
		sumvy = sumvy + vy[i];
		sumvz = sumvz + vz[i];
	}

	sumvx = sumvx / (float)(n);
	sumvy = sumvy / (float)(n);
	sumvz = sumvz / (float)(n);

	for (i = 0; i < n; i++)
	{
		vx[i] = (vx[i] - sumvx) * 1e-5;
		vy[i] = (vy[i] - sumvy) * 1e-5;
		vz[i] = (vz[i] - sumvz) * 1e-5;
	}

	// total edge of the cube
	float boxlx = a * NC;
	float boxly = a * NC;
	float boxlz = a * NC;

	//center the cube in (0,0,0)
	for (i = 0; i < n; i++)
	{
		rx[i] = rx[i] - boxlx * 0.5;
		ry[i] = ry[i] - boxly * 0.5;
		rz[i] = rz[i] - boxlz * 0.5;
	}

	//write the positions in file
	for (i = 0; i < n; i++)
		fprintf(in, "12 %e %e %e %e %e %e \n", rx[i], ry[i], rz[i], vx[i], vy[i], vz[i]);
	fclose(in);
	return 0;
}

//determine a random seed for generate the initial velocities	
float ranf(int seed[1])
{
	int l = 1029;
	int c = 221591;
	int m = 1048576;

	seed[0] = (seed[0] * l + c) % m;
	float fun = (1.0 * seed[0]) / (1.0 * m);
	return fun;
}

float gauss(int seed[1])
{
	int i;
	float a1 = 3.949846138;
	float a3 = 0.252408784;
	float a5 = 0.076542912;
	float a7 = 0.008355968;
	float a9 = 0.029899776;
	float sum = 0.0;
	for (i = 0; i < 12; i++)
		sum = sum + ranf(seed);

	float r = 0.25 * (sum - 6.0);
	float r2 = r * r;
	float fun = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r;
	return fun;
}
