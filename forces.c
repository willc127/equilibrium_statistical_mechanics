#include <stdio.h>
#include <math.h>
#define max 10000

void accel(int n, float eps24, float rcutsq, float rx[max], float ry[max], float rz[max], float ax[max], float ay[max], float az[max], float sigsq, float boxlx, float boxly, float boxlz);
float potential(float rijsq, float sigsq, float *wij);

int main()
{
	int i, n;
	float rx[max], ry[max], rz[max], vx[max], vy[max], vz[max], ax[max], ay[max], az[max];
	char dummy; // to read the empty line in cfc file
	FILE *in;

	int atom = 12;
	n = 500; // number of particles
	in = fopen("cfc.xyz", "r");
	fscanf(in, "%d", &n);
	printf("%d\n", n);
	fscanf(in, "%c", &dummy);

	for (i = 0; i < n; i++)
	{
		fscanf(in, "%d %f %f %f %f %f %f", &atom, &rx[i], &ry[i], &rz[i], &vx[i], &vy[i], &vz[i]);
	}

	fclose(in);
	float m = 16.0; // molecular weight
	float NA = 6.0221409E23;
	float dens = 100.0;												 // density in kg/m3
	float a = pow(4.0 * m / (1000.0 * dens * NA), 1.0 / 3.0) * 1e10; // edge em angstrom
	float kb = 1.38064852E-23;
	float r = NA * kb;
	float sigma = 3.73; // angstrom
	float eps = 148.0;	// in Kelvin
	float eps24 = 24 * eps * r / (m / 1000.0);
	float rcutsq = 2.5 * sigma * 2.5 * sigma;
	int NC = cbrt(n * 0.25); // number of cells
	printf("nc = %d\n", NC);
	float sigsq = sigma * sigma;

	// define the total edge length of cube
	float boxlx = a * NC;
	float boxly = a * NC;
	float boxlz = a * NC;

	for (i = 0; i < n; i++)
	{ // centralizar o cubo na posi��o (0,0,0)
		rx[i] = rx[i] - boxlx * 0.5;
		ry[i] = ry[i] - boxly * 0.5;
		rz[i] = rz[i] - boxlz * 0.5;
	}

	accel(n, eps24, rcutsq, rx, ry, rz, ax, ay, az, sigsq, boxlx, boxly, boxlz);

	FILE *out;
	out = fopen("acel.xyz", "w");
	printf("n = %d\n", n);
	fprintf(in, "%d\n\n", n);
	for (i = 0; i < n; i++)
		fprintf(out, "12 %e %e %e %e %e %e %e %e %e \n", rx[i], ry[i], rz[i], vx[i], vy[i], vz[i], ax[i], ay[i], az[i]);
	fclose(out);

	return 0;
}

void accel(int n, float eps24, float rcutsq, float rx[max], float ry[max], float rz[max], float ax[max], float ay[max], float az[max], float sigsq, float boxlx, float boxly, float boxlz)
{

	float rxi, ryi, rzi, wij, aij, rxij, ryij, rzij, pot, rijsq;
	int i, j;

	for (i = 0; i < n; i++)
	{
		ax[i] = 0.0;
		ay[i] = 0.0;
		az[i] = 0.0;
	}

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
				pot = potential(rijsq, sigsq, &wij);
				aij = 1.0 * (eps24 * (wij / rijsq)) * 1e-10; // A/fs�
				ax[i] = ax[i] + aij * rxij;
				ay[i] = ay[i] + aij * ryij;
				az[i] = az[i] + aij * rzij;
				ax[j] = ax[j] - aij * rxij;
				ay[j] = ay[j] - aij * ryij;
				az[j] = az[j] - aij * rzij;
			}
		}
	}
}

float potential(float rijsq, float sigsq, float *wij)
{
	float sr2 = sigsq / rijsq;
	float sr6 = sr2 * sr2 * sr2;
	float sr12 = sr6 * sr6;
	float pot = sr12 - sr6;
	*wij = 2.0 * sr12 - sr6;
	return pot;
}
