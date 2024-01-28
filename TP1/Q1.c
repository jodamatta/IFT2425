#include <stdio.h>
#include <math.h>

#define EPSILON 1e-5
#define MAX_ITERATIONS 200

double f(double y[], int N, double c) {
	double numerator = 0.0;
	double denominator = 0.0;
	double sum_log_y = 0.0;

	for (int i = 0; i < N; ++i) {
		numerator += pow(y[i], c) * log(y[i]);
		denominator += pow(y[i], c);
		sum_log_y += log(y[i]);
	}
	return numerator / denominator - 1.0 / c - sum_log_y / N;
}

double derivativeA(double y[], int N, double c, double epsilon) {
	return (f(y, N, c + epsilon) - f(y, N, c)) / epsilon;
}

double methode_newtonA(double y[], int N, double initial_guess, double tol) {
	double c_n = initial_guess;

	for (int i = 0; i < 100; ++i) {
		double f_c_n = f(y, N, c_n);
		double f_prime_c_n = derivativeA(y, N, c_n, EPSILON);
		c_n = c_n - f_c_n / f_prime_c_n;

		if (fabs(f_c_n) < tol) {
			break;
		}
	}
	return c_n;
}

double derivativeB(double y[], int N, double c, double epsilon) {
	return (-f(y, N, c + 2 * epsilon) + 8 * f(y, N, c + epsilon) - 8 * f(y, N, c - epsilon) + f(y, N, c - 2 * epsilon)) / (12 * epsilon);
}

double methode_newtonB(double y[], int N, double initial_guess, double tol) {
	double c_n = initial_guess;

	for (int i = 0; i < 100; ++i) {
		double f_c_n = f(y, N, c_n);
		double f_prime_c_n = derivativeB(y, N, c_n, EPSILON);
		c_n = c_n - f_c_n / f_prime_c_n;

		if (fabs(f_c_n) < tol) {
			break;
		}
	}
	return c_n;
}

double derivee_analytique(double y[], int N, double c) {
	// trouve avec wolfram alpha
	return 0.254693 -  1 / c;
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
	
	double root = methode_newtonA(y_values, N, initial_guess, tol);
	printf("Q1A: racine trouve: %lf\n", root);

	root = methode_newtonB(y_values, N, initial_guess, tol);
	printf("Q1B: racine trouve: %lf\n", root);

	root = methode_newtonC(y_values, N, initial_guess, tol);
	printf("Q1C: racine trouve: %lf\n", root);
}