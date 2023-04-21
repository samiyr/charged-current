#ifndef TESTS_H
#define TESTS_H

#include "DIS.cpp"
#include <iostream>
#include "Utility.cpp"

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
	PDFInterface pdf("CT18NLO", 0);

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

void pdf_comparison_tests() {
	PDFInterface pdf1("CT18NLO", 0);
	PDFInterface pdf2("CT18NLO", 0);

	const double x = 0.1;
	const double Q2 = 200;
	const double tol = 1e-15;
	pdf1.evaluate(x, Q2);

	double_comparison(pdf1.xf(Flavor::Up), pdf2.xf_evaluate(Flavor::Up, x, Q2), tol);
	double_comparison(pdf1.xf(Flavor::Down), pdf2.xf_evaluate(Flavor::Down, x, Q2), tol);
	double_comparison(pdf1.xf(Flavor::Charm), pdf2.xf_evaluate(Flavor::Charm, x, Q2), tol);
	double_comparison(pdf1.xf(Flavor::Strange), pdf2.xf_evaluate(Flavor::Strange, x, Q2), tol);
	double_comparison(pdf1.xf(Flavor::Top), pdf2.xf_evaluate(Flavor::Top, x, Q2), tol);
	double_comparison(pdf1.xf(Flavor::Bottom), pdf2.xf_evaluate(Flavor::Bottom, x, Q2), tol);
	double_comparison(pdf1.xg(), pdf2.xf_evaluate(Flavor::Gluon, x, Q2), tol);
}

void ckm_tests() {
	const double tol = 1e-12;

	double_comparison(CKM::squared(Flavor::Up, Flavor::Down), Constants::V_ud, tol);
	double_comparison(CKM::squared(Flavor::Down, Flavor::Up), Constants::V_ud, tol);
	
	double_comparison(CKM::squared(Flavor::Up, Flavor::Strange), Constants::V_us, tol);
	double_comparison(CKM::squared(Flavor::Strange, Flavor::Up), Constants::V_us, tol);
	
	double_comparison(CKM::squared(Flavor::Up, Flavor::Bottom), Constants::V_ub, tol);
	double_comparison(CKM::squared(Flavor::Bottom, Flavor::Up), Constants::V_ub, tol);


	double_comparison(CKM::squared(Flavor::Charm, Flavor::Down), Constants::V_cd, tol);
	double_comparison(CKM::squared(Flavor::Down, Flavor::Charm), Constants::V_cd, tol);
	
	double_comparison(CKM::squared(Flavor::Charm, Flavor::Strange), Constants::V_cs, tol);
	double_comparison(CKM::squared(Flavor::Strange, Flavor::Charm), Constants::V_cs, tol);
	
	double_comparison(CKM::squared(Flavor::Charm, Flavor::Bottom), Constants::V_cb, tol);
	double_comparison(CKM::squared(Flavor::Bottom, Flavor::Charm), Constants::V_cb, tol);


	double_comparison(CKM::squared(Flavor::Top, Flavor::Down), Constants::V_td, tol);
	double_comparison(CKM::squared(Flavor::Down, Flavor::Top), Constants::V_td, tol);
	
	double_comparison(CKM::squared(Flavor::Top, Flavor::Strange), Constants::V_ts, tol);
	double_comparison(CKM::squared(Flavor::Strange, Flavor::Top), Constants::V_ts, tol);
	
	double_comparison(CKM::squared(Flavor::Top, Flavor::Bottom), Constants::V_tb, tol);
	double_comparison(CKM::squared(Flavor::Bottom, Flavor::Top), Constants::V_tb, tol);

	const std::vector<FlavorType> upper = {Flavor::Up, Flavor::Charm, Flavor::Top};
	for (int i = 0; i < upper.size(); i++) {
		for (int j = 0; j < upper.size(); j++) {
			double_comparison(CKM::squared(upper[i], upper[j]), 0, tol);
		}
	}
	const std::vector<FlavorType> lower = {Flavor::Down, Flavor::Strange, Flavor::Bottom};
	for (int i = 0; i < lower.size(); i++) {
		for (int j = 0; j < lower.size(); j++) {
			double_comparison(CKM::squared(lower[i], lower[j]), 0, tol);
		}
	}
}

int main(int argc, char const *argv[]) {
	integration_tests();
	coefficient_tests();
	pdf_evaluation_tests();
	
	return 0;
}


#endif