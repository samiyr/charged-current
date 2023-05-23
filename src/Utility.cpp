#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>

#define POW2(x) (x) * (x)
#define POW4(x) (x) * (x) * (x) * (x)

template <typename T>
constexpr std::vector<T> vector_intersection(std::vector<T> &v1, std::vector<T> &v2) {
    std::vector<T> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(v3));
    return v3;
}

bool double_comparison(double a, double b, double tolerance = 1e-5) {
	bool flag = abs(a - b) < tolerance;
	// const bool flag = std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;

	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	if (flag) {
		std::cout << "PASS" << std::endl;
	} else {
		std::cout << "FAIL" << std::endl;
	}
	return flag;
}
bool double_comparison_rel(double a, double b, double tolerance = 1e-5) {
	const bool flag = std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;

	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	if (flag) {
		std::cout << "PASS" << std::endl;
	} else {
		std::cout << "FAIL" << std::endl;
	}
	return flag;
}

double order_of_magnitude(const double x) {
	return std::log10(std::abs(x));
}

bool double_comparison_digits(double a, double b, int digits = 5) {
	const double difference = order_of_magnitude(a - b);
	const double original = order_of_magnitude(a);
	bool flag = abs(difference - original) >= double(digits);

	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	if (flag) {
		std::cout << "PASS" << std::endl;
	} else {
		std::cout << "FAIL" << std::endl;
	}
	return flag;
}


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace Utility {
	constexpr double logm1(const double x) {
		return std::log1p(-x);
	}

	constexpr double gamma(const double z) {
		return gsl_sf_gamma(z);
	}

	constexpr double beta(const double a, const double b) {
		return boost::math::beta(a, b);
	}
	constexpr double incomplete_beta(const double z, const double a, const double b) {
		return boost::math::beta(a, b, z);
	}
	// inline double beta(const double a, const double b) {
	// 	return gsl_sf_beta(a, b);
	// }
	// inline double incomplete_beta(const double z, const double a, const double b) {
	// 	return gsl_sf_beta_inc(a, b, z) * gsl_sf_beta(a, b);
	// }

	void non_aborting_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
		std::cout << "gsl error " << gsl_errno << " (" << std::string(file) << "): " << std::string(reason) << std::endl;
	}
}

#endif