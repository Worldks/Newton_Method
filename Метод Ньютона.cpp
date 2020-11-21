// 10 Вариант
#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

double const E = 0.000000001;

void Gaus(double* const A[], double const B[], const int n, double x[]);

double func1(double x[]);

double func2(double x[]);

double dx1_func1(double x[]);

double dx2_func1(double x[]);

double dx1_func2(double x[]);

double dx2_func2(double x[]);

void Jacobi(double** J, double x[], double M);

void Newton(double x[], int n, int NIT, double M = 0);



int main()
{
	int const n = 2;
	double x[n];
	x[0] = 1;
	x[1] = 1;
	int NIT = 100; // предельное число итераций
	//(1;1)
	cout << "(1;1):" << endl;
	cout << "Bez M:\n";
	Newton(x, n, NIT);

	x[0] = 1;
	x[1] = 1;
	cout << "M = 0.01:\n";
	Newton(x, n, NIT, 0.01);

	x[0] = 1;
	x[1] = 1;
	cout << "M = 0.05:\n";
	Newton(x, n, NIT, 0.05);

	x[0] = 1;
	x[1] = 1;
	cout << "M = 0.1:\n";
	Newton(x, n, NIT, 0.1);

	system("pause");
	return 0;
}

void Gaus(double* const A[]/*Якобиан*/, double const B[]/*Вектор невязки*/, const int n, double x[]) {

	double** a = new double* [n];
	for (int i = 0; i < n; i++) {
		a[i] = new double[n];
		for (int j = 0; j < n; j++)
			a[i][j] = A[i][j];
	}

	double* b = new double[n];
	for (int i = 0; i < n; i++)
		b[i] = B[i];

	for (int k = 0; k < n - 1; k++) {
		double AKK = 0;
		int p = k;

		for (int i = k; i < n; i++) {
			if (fabs(a[i][k]) < AKK)
				continue;
			p = i;
			AKK = fabs(a[i][k]);
		}


		if (p != k) {
			double* c = a[p];
			a[p] = a[k];
			a[k] = c;
			double d = b[p];
			b[p] = b[k];
			b[k] = d;
		}

		double AMAIN = a[k][k];
		if (!AMAIN)
			return;

		for (int j = k; j < n; j++)
			a[k][j] = a[k][j] / AMAIN;
		b[k] = b[k] / AMAIN;
		for (int i = k + 1; i < n; i++) {
			for (int j = k + 1; j < n; j++)
				a[i][j] = a[i][j] - a[i][k] * a[k][j];
			b[i] = b[i] - a[i][k] * b[k];
		}

		/*for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				cout << setw(10) << a[i][j];
			cout << "\n\n";
		}*/
	}

	if (!a[n - 1][n - 1])
		return;
	b[n - 1] = b[n - 1] / a[n - 1][n - 1];
	x[n - 1] = b[n - 1];

	for (int k = n - 2; k >= 0; k--) {
		double summ = 0;
		for (int j = k + 1; j < n; j++)
			summ += a[k][j] * x[j];
		x[k] = b[k] - summ;
	}

	/*cout << "X:\n";
	for (int i = 0; i < n; i++)
		cout << x[i] << "\n";
	*/

	/*double* F = new double[n];
	//A * x - b
	for (int i = 0; i < n; i++) {
		F[i] = 0;
		for (int j = 0; j < n; j++)
			F[i] += A[i][j] * x[j];
		F[i] -= B[i];
	}
	double maxF = F[0];
	cout << "F:\n";
	for (int i = 0; i < n; i++)
		cout << F[i] << "\n";*/
}

double func1(double x[]) {
	return (sin(x[0] + 1) - x[1] - 1);
}

double func2(double x[]) {
	return (2*x[0] + cos(x[1]) - 2);
}

double dx1_func1(double x[]) {
	return cos(x[0]);
}

double dx2_func1(double x[]) {
	return -1;
}

double dx1_func2(double x[]) {
	return 2;
}

double dx2_func2(double x[]) {
	return -1*sin(x[1]);
}


void Jacobi(double** J, double x[], double M = 0) {
	if (M == 0) {
		J[0][0] = dx1_func1(x); J[0][1] = dx2_func1(x);
		J[1][0] = dx1_func2(x); J[1][1] = dx2_func2(x);
	}
	else {  // Вычисление матрицы Якоби конечно-разностным методом
		double x1[2];
		x1[0] = x[0];
		x1[1] = x[1];

		x1[0] += x[0] * M;
		J[0][0] = (func1(x1) - func1(x)) / (x1[0] - x[0]);

		x1[0] = x[0];
		x1[1] = x[1];

		x1[1] += x[1] * M;
		J[0][1] = (func1(x1) - func1(x)) / (x[1] * M);

		x1[0] = x[0];
		x1[1] = x[1];

		x1[0] += x[0] * M;
		J[1][0] = (func2(x1) - func2(x)) / (x[0] * M);

		x1[0] = x[0];
		x1[1] = x[1];

		x1[1] += x[1] * M;
		J[1][1] = (func2(x1) - func2(x)) / (x[1] * M);
	}
}

void Newton(double x[], int n, int NIT, double M) {
	int k = 1;
	double d1, d2;   // Необходимы для определения условия выхода из цикла e1 and e2 ТОЕСТЬ  E, т.к. они равны

	double* dx = new double[n];
	double* F = new double[n];
	double** J = new double* [n];
	for (int i = 0; i < n; i++)
		J[i] = new double[n];

	while (true) {
		cout << "#" << k;
		if (k > 1)
			cout << "    " << setw(15) << d1 << setw(15) << d2;
		cout << "\n";

		F[0] = -func1(x);
		F[1] = -func2(x);

		Jacobi(J, x, M); // После этой функции изменяется Матрица Якобиана

		Gaus(J, F, n, dx);  // После этой функции изменяется вектор x и dx являющийся вектором поправки

		double* x_old = new double[n]; //Массив сохраняющий значение вектора неизвестных до уточнения решения(приближение)
		for (int i = 0; i < n; i++) {
			x_old[i] = x[i];
			x[i] += dx[i];  // Уточнение решения
		}

		d1 = fabs(func1(x));
		if (d1 < fabs(func2(x)))
			d1 = fabs(func2(x));

		if (fabs(x[0]) < 1)
			d2 = fabs(x[0] - x_old[0]);
		else
			d2 = fabs((x[0] - x_old[0]) / x[0]);

		if (fabs(x[1]) < 1) {
			if (d2 < fabs(x[1] - x_old[1]))
				d2 = fabs(x[1] - x_old[1]);
		}
		else
			if (d2 < fabs((x[1] - x_old[1]) / x[1]))
				d2 = fabs((x[1] - x_old[1]) / x[1]);

		if (d1 < E && d2 < E)
			break;

		if (k >= NIT) {
			cout << "IER = 2";
			return;
		}

		k++;
	}

	cout << "X: " << x[0] << "   " << x[1] << "\n\n\n";
}