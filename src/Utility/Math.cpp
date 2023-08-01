#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <numeric>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>

namespace Math {
	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}
	
	constexpr double log1m(const double x) {
		return std::log(1 - x);
		// return std::log1p(-x);
	}

	inline double gamma(const double z) {
		return gsl_sf_gamma(z);
	}

	inline double beta(const double a, const double b) {
		return boost::math::beta(a, b);
	}
	inline double incomplete_beta(const double z, const double a, const double b) {
		// return std::pow(z, a) * (1.0 / a + z * (1 - b) / (1 + a) + z * z * (2 - b) * (1 - b) / (2 * (2 + a)));
		return boost::math::beta(a, b, z);
	}
	// inline double beta(const double a, const double b) {
	// 	return gsl_sf_beta(a, b);
	// }
	// inline double incomplete_beta(const double z, const double a, const double b) {
	// 	return gsl_sf_beta_inc(a, b, z) * gsl_sf_beta(a, b);
	// }

	template <typename Iterator1, typename Iterator2>
	constexpr double dot_product(const Iterator1 &iter1, const Iterator2 &iter2) {
		return std::inner_product(std::begin(iter1), std::end(iter1), std::begin(iter2), 0.0);
	}

	constexpr double kahan_sum(const std::initializer_list<double> values) {
		double sum = 0.0;
		double c = 0.0;
		for (const double value : values) {
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}
		return sum;
	}
}

#endif