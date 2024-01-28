#include <stdio.h>
#include <math.h>

#define EPSILON 1e-5

double derivee_analytique(double y[], int N, double c) {
	// trouve avec wolfram alpha
	return 0.254693 - 1 / c;
}

double methode_newtonC(double y[], int N, double c, double tol) {
	double c_n = c;
	double value;

	for (int i = 0; i < 100; ++i) {
		double f_c_n = f(y, N, c_n);
		double f_prime_c_n = derivee_analytique(y, N, c_n);
		c_n = c_n - f_c_n / f_prime_c_n;

		if (fabs(f_c_n) < tol) {
			break;
		}
	}
	return c_n;
}

int main() {

	int N = 10;
	double y_values[] = { 0.11,0.24,0.27,0.52,1.13,1.54,1.71,1.84,1.92,2.01 };

	double initial_guess = 0.25;
	double tol = 1e-6;

	double root = methode_newtonC(y_values, N, initial_guess, tol);
	printf("Q1C: racine trouve: %lf\n", root);
}