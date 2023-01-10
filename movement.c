#include <stdio.h>
#include <math.h>
#define max 10000

void accel(int n, float eps24, float rcutsq, float rx[max], float ry[max], float rz[max], float ax[max], float ay[max], float az[max], float sigsq, float boxlx, float boxly, float boxlz);
float potential(float rijsq, float sigsq, float *wij);

int main()
{
	int i, n, k;
	float rx[max], ry[max], rz[max], vx[max], vy[max], vz[max], ax[max], ay[max], az[max];
	char dummy; // para ler a linha vazia do arquivo cfc
	FILE *in;

	int atom = 12;				// part�cula a ser utilizada
	n = 500;					// n�mero de part�culas
	in = fopen("cfc.xyz", "r"); // usar'r' somente para leitura
	fscanf(in, "%d", &n);
	printf("%d\n", n);
	fscanf(in, "%c", &dummy);

	for (i = 0; i < n; i++)
	{
		fscanf(in, "%d %f %f %f %f %f %f", &atom, &rx[i], &ry[i], &rz[i], &vx[i], &vy[i], &vz[i]); // o caractere & serve para salvar o valor da vari�vel a cada rodada
	}
	fclose(in);

	float m = 16.0; // massa molar
	float NA = 6.0221409E23;
	float dens = 100.0;												 // Densidade em kg/m�
	float a = pow(4.0 * m / (1000.0 * dens * NA), 1.0 / 3.0) * 1e10; // aresta em angstrom

	float kb = 1.38064852E-23;
	float r = NA * kb;
	float sigma = 3.73; // angstrom
	float eps = 148.0;	// em K
	float eps24 = 24 * eps * r / (m * 1e-3);
	float rcutsq = 2.5 * sigma * 2.5 * sigma;
	int NC = cbrt(n * 0.25); // n�mero de c�lulas
	printf("nc = %d\n", NC);
	float sigsq = sigma * sigma;

	float boxlx = a * NC; // define a aresta total do cubo
	float boxly = a * NC;
	float boxlz = a * NC;

	accel(n, eps24, rcutsq, rx, ry, rz, ax, ay, az, sigsq, boxlx, boxly, boxlz);

	FILE *out;
	out = fopen("movs.xyz", "w");

	float t = 0;
	float dt = 1;
	float dt2 = dt * 0.5;
	float tf = 5000;
	float T;
	float kin;
	int ct = 0;

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

		accel(n, eps24, rcutsq, rx, ry, rz, ax, ay, az, sigsq, boxlx, boxly, boxlz);
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
	}

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
				aij = wij / rijsq; // A/fs�
				ax[i] = ax[i] + aij * rxij;
				ay[i] = ay[i] + aij * ryij;
				az[i] = az[i] + aij * rzij;
				ax[j] = ax[j] - aij * rxij;
				ay[j] = ay[j] - aij * ryij;
				az[j] = az[j] - aij * rzij;
			}
		}
	}

	for (i = 0; i < n; i++)
	{
		ax[i] = ax[i] * eps24 * 1e-10;
		ay[i] = ay[i] * eps24 * 1e-10;
		az[i] = az[i] * eps24 * 1e-10;
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
