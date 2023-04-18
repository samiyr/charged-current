#ifndef TESTS_H
#define TESTS_H

#include "DIS.cpp"
#include <iostream>

void double_comparison(double a, double b, double tolerance = 1e-5) {
	bool flag = abs(a - b) < tolerance;

	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	if (flag) {
		std::cout << "PASS" << std::endl;
	} else {
		std::cout << "FAIL" << std::endl;
	}
}

double test_integrand_1(double x[], size_t dim, void *params) {
	return 1.0 / std::sqrt(x[0]);
}

double test_integrand_2(double x[], size_t dim, void *params) {
	return 1.0 / std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

void coefficient_tests() {
	std::cout << "DIS C2q test" << std::endl;
	double_comparison(Coefficients::DIS::C2q(0.001), 6.31063, 1e-5);
	double_comparison(Coefficients::DIS::C2q(0.01), 4.61014, 1e-5);
	double_comparison(Coefficients::DIS::C2q(0.1), 2.99933, 1e-5);
	double_comparison(Coefficients::DIS::C2q(0.5), 0.75, 1e-5);
	double_comparison(Coefficients::DIS::C2q(0.99), -791.17476, 1e-5);
	std::cout << std::endl << std::endl;

	std::cout << "DIS C2g test" << std::endl;
	double_comparison(Coefficients::DIS::C2g(0.001), 2.95047, 1e-5);
	double_comparison(Coefficients::DIS::C2g(0.01), 1.79167, 1e-5);
	double_comparison(Coefficients::DIS::C2g(0.1), 0.760862, 1e-5);
	double_comparison(Coefficients::DIS::C2g(0.5), 0.5, 1e-5);
	double_comparison(Coefficients::DIS::C2g(0.99), -2.71247, 1e-5);
	std::cout << std::endl << std::endl;
}

void integration_tests() {
	Integrator int1(&test_integrand_1, {0.001}, {1}, 10'000);
	Integrator::Result res1 = int1.integrate();
	
	std::cout << "Integration test #1" << std::endl;
	std::cout << res1 << std::endl;
	double_comparison(res1.value, 1.93675, 1e-5);
	std::cout << std::endl << std::endl;
	
	Integrator int2(&test_integrand_2, {0.001, 0.001, 0.001}, {1, 1, 1}, 200'000);
	Integrator::Result res2 = int2.integrate();
	
	std::cout << "Integration test #2" << std::endl;
	std::cout << res2 << std::endl;
	double_comparison(res2.value, 1.18478, 1e-5);
	std::cout << std::endl << std::endl;
}

void pdf_evaluation_tests() {
	PDFEvaluation pdf("CT18NLO", 0);

	std::cout << "CT18NLO central PDF test with up at Q^2 = 4 GeV^2:" << std::endl;
	double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-3, 4), 0.488643);
	double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-2, 4), 0.444769);
	double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-1, 4), 0.64825);
	double_comparison(pdf.xf_evaluate(Flavor::Up, 0.2, 4), 0.700196);
	double_comparison(pdf.xf_evaluate(Flavor::Up, 0.5, 4), 0.287891);
	double_comparison(pdf.xf_evaluate(Flavor::Up, 0.9, 4), 0.00186158);
	
	std::cout << "CT18NLO central PDF test with down at Q^2 = 4 GeV^2:" << std::endl;

	double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-3, 4), 0.447382);
	double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-2, 4), 0.395838);
	double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-1, 4), 0.438871);
	double_comparison(pdf.xf_evaluate(Flavor::Down, 0.2, 4), 0.369744);
	double_comparison(pdf.xf_evaluate(Flavor::Down, 0.5, 4), 0.0831107);
	double_comparison(pdf.xf_evaluate(Flavor::Down, 0.9, 4), 0.000204438);
}

int main(int argc, char const *argv[]) {
	integration_tests();
	coefficient_tests();
	pdf_evaluation_tests();
	
	return 0;
}


#endif