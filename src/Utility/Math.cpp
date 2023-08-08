#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <numeric>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>
#include <concepts>

namespace Math {
	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	constexpr double log1m(const double x) {
		return std::log(1 - x);
	}

	inline double gamma(const double z) {
		return gsl_sf_gamma(z);
	}

	inline double beta(const double a, const double b) {
		return boost::math::beta(a, b);
	}
	inline double incomplete_beta(const double z, const double a, const double b) {
		return boost::math::beta(a, b, z);
	}

	template <typename Iterator1, typename Iterator2>
	constexpr double dot_product(const Iterator1 &iter1, const Iterator2 &iter2) {
		return std::inner_product(std::begin(iter1), std::end(iter1), std::begin(iter2), 0.0);
	}

	template <typename T>
	concept NumericType= std::integral<T> || std::floating_point<T>;

	template <NumericType T>
	static constexpr T pow2(const T x) {
		return x * x;
	}

	template <NumericType T>
	static constexpr T pow4(const T x) {
		return x * x * x * x;
	}
}

#endif