#include <stdio.h>
#include <math.h>
#define max 10000

void accel(float m, float kb, int n, float eps24, float eps, float rcutsq,
		   float rx[max], float ry[max], float rz[max], float ax[max], float ay[max], float az[max],
		   float sigsq, float boxlx, float boxly, float boxlz, float *V, float *W);

float potential(float rijsq, float sigsq, float *wij, float eps24, float kb);

int main()
{
	int i, n, k;
	float rx[max], ry[max], rz[max], vx[max], vy[max], vz[max], ax[max], ay[max], az[max];
	char dummy;
	FILE *in;

	int atom = 12; // define particle
	n = 864;	   // total of particles

	// read the initial configuration
	in = fopen("cfc.xyz", "r");
	fscanf(in, "%d", &n);
	printf("%d\n", n);
	fscanf(in, "%c", &dummy); // read the empty row of cfc.xyz file

	for (i = 0; i < n; i++)
	{
		fscanf(in, "%d %f %f %f %f %f %f", &atom, &rx[i], &ry[i], &rz[i], &vx[i], &vy[i], &vz[i]);
	}
	fclose(in);

	float m = 16.0;													 // molecular weigth (methane)
	float NA = 6.0221409E23;										 // Avogadro's number
	float dens = 100.0;												 // Density in kg/m3
	float a = pow(4.0 * m / (1000.0 * dens * NA), 1.0 / 3.0) * 1e10; // edge in angstrom
	float kb = 1.38064852E-23;										 // Boltzmann's constant
	float r = NA * kb;												 // universal gas constant
	float sigma = 3.73;												 // angstrom
	float eps = 148.0;												 // in Kelvin
	float eps24 = 24 * eps * r / (m * 1e-3);
	float rcutsq = 2.5 * sigma * 2.5 * sigma;
	int NC = cbrt(n * 0.25); // number of cells
	float sigsq = sigma * sigma;
	float boxlx = a * NC; // total edge of the cube x-direction
	float boxly = a * NC; // total edge of the cube y-direction
	float boxlz = a * NC; // total edge of the cube z-direction
	float V;
	float W;
	float wij;

	accel(m, kb, n, eps24, eps, rcutsq, rx, ry, rz, ax, ay, az, sigsq, boxlx, boxly, boxlz, &V, &W);

	// open file to store the position of each atom over time
	FILE *out;
	out = fopen("movs.xyz", "w");

	// Define time range
	float t = 0;		  // initial time
	float dt = 1;		  // step
	float dt2 = dt * 0.5; // half-step
	float tf = 3000;	  // final time

	float T;
	float kin;
	float pot;

	// open file to store the thermodynamical variables
	FILE *temp;
	temp = fopen("temper.dat", "w");

	while (t < tf)
	{
		printf("Time = %f\n", t);
		for (i = 0; i < n; i++)
		{
			vx[i] = vx[i] + ax[i] * dt2;
			vy[i] = vy[i] + ay[i] * dt2;
			vz[i] = vz[i] + az[i] * dt2;

			rx[i] = rx[i] + vx[i] * dt;
			ry[i] = ry[i] + vy[i] * dt;
			rz[i] = rz[i] + vz[i] * dt;
		}

		accel(m, kb, n, eps24, eps, rcutsq, rx, ry, rz, ax, ay, az, sigsq, boxlx, boxly, boxlz, &V, &W);

		float sumv2 = 0.0;

		for (i = 0; i < n; i++)
		{
			vx[i] = vx[i] + ax[i] * dt2;
			vy[i] = vy[i] + ay[i] * dt2;
			vz[i] = vz[i] + az[i] * dt2;
			sumv2 = sumv2 + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
		}
		fprintf(out, "%d\n\n", n);

		for (k = 0; k < n; k++)
		{
			fprintf(out, "12 %f %f %f  \n", rx[k], ry[k], rz[k]);
		}

		t = t + dt;
		float rijsq;
		float sigsq;
		float kin = 0.5 * m * sumv2 * 1e7 / (float)(n); // J/molecule
		float H = kin + V;
		float T = 2.0 * kin / (3.0 * r);
		float Vol = a * NC * a * NC * a * NC * 1e-30;											 // m3
		float Pr = ((float)(n) / NA * r * T / Vol + W / 3 * eps24 * m * 1e-3 / Vol / NA) * 1e-5; // bar
		fprintf(temp, "%f %f %f %f %f %f\n", t, T, Pr, kin, V, H);
	}

	fclose(out);
	return 0;
}

// Create a function to calculate the acceleration vector
void accel(float m, float kb, int n, float eps24, float eps, float rcutsq,
		   float rx[max], float ry[max], float rz[max], float ax[max], float ay[max], float az[max],
		   float sigsq, float boxlx, float boxly, float boxlz, float *V, float *W)
{

	float rxi, ryi, rzi, aij, rxij, ryij, rzij, pot, rijsq, wij;
	int i, j;

	for (i = 0; i < n; i++)
	{
		ax[i] = 0.0;
		ay[i] = 0.0;
		az[i] = 0.0;
	}

	*V = 0.0;
	*W = 0.0;
	for (i = 0; i < n - 1; i++)
	{
		rxi = rx[i];
		ryi = ry[i];
		rzi = rz[i];
		for (j = i + 1; j < n; j++)
		{
			rxij = rxi - rx[j];
			ryij = ryi - ry[j];
			rzij = rzi - rz[j];
			rxij = rxij - boxlx * round(rxij / boxlx);
			ryij = ryij - boxly * round(ryij / boxly);
			rzij = rzij - boxlz * round(rzij / boxlz);
			rijsq = rxij * rxij + ryij * ryij + rzij * rzij;

			if (rijsq <= rcutsq)
			{
				pot = potential(rijsq, sigsq, &wij, eps, kb);
				*V = *V + pot;
				aij = wij / rijsq * 1e-10; // A/fs�
				ax[i] = ax[i] + aij * rxij;
				ay[i] = ay[i] + aij * ryij;
				az[i] = az[i] + aij * rzij;
				ax[j] = ax[j] - aij * rxij;
				ay[j] = ay[j] - aij * ryij;
				az[j] = az[j] - aij * rzij;
				*W = *W + wij;
			}
		}
	}
	*V = *V * eps24 * m * 1e-3 / 6.0 / (float)(n);
	for (i = 0; i < n; i++)
	{
		ax[i] = ax[i] * eps24;
		ay[i] = ay[i] * eps24;
		az[i] = az[i] * eps24;
	}
}

// Create a function to calculate the potential energy
float potential(float rijsq, float sigsq, float *wij, float eps, float kb)
{
	float sr2 = sigsq / rijsq;
	float sr6 = sr2 * sr2 * sr2;
	float sr12 = sr6 * sr6;
	float pot = sr12 - sr6;

	*wij = 2.0 * sr12 - sr6;

	return pot;
}
