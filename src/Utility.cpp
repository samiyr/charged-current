#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <numeric>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>
#include <string>
#include <thread>
#include <limits>

#define POW2(x) (x) * (x)
#define POW4(x) (x) * (x) * (x) * (x)

namespace IO {
	inline std::ostream& endl(std::ostream& os) {
		os.put(os.widen('\n'));
		return os;
	}
}

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
		std::cout << "PASS" << IO::endl;
	} else {
		std::cout << "FAIL" << IO::endl;
	}
	return flag;
}
bool double_comparison_rel(double a, double b, double tolerance = 1e-5) {
	const bool flag = std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;

	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	if (flag) {
		std::cout << "PASS" << IO::endl;
	} else {
		std::cout << "FAIL" << IO::endl;
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
		std::cout << "PASS" << IO::endl;
	} else {
		std::cout << "FAIL" << IO::endl;
	}
	return flag;
}


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace Utility {
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

	void non_aborting_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
		std::cout << "gsl error " << gsl_errno << " (" << std::string(file) << ", " << line << "): " << std::string(reason) << IO::endl;
	}

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

	constexpr std::string bool_to_string(const bool input) {
		return input ? "true" : "false";
	}

	// https://stackoverflow.com/a/19743821/1531270
	template<class Class>
	struct Traced {
	public:
		Traced() = default;
		Traced(Traced const&) {
			#pragma omp critical
			std::cout << typeid(Class).name() << " copy constructor called" << IO::endl;
		}
	protected:
		~Traced() = default;
	};

	static unsigned int get_default_thread_count() {
		const unsigned int hardware = std::thread::hardware_concurrency();
		return std::max(1U, hardware);
	}

	template <typename T>
	static void multiply(std::vector<T> &v1, const std::vector<T> &v2) {
		for (typename std::vector<T>::size_type i = 0; i < v1.size(); i++) {
			v1[i] *= v2[i];
		}
	}

	constexpr int size_to_int(const size_t input) {
		if (input > std::numeric_limits<int>::max()) {
			throw std::runtime_error("Cannot cast size_t larger than std::numeric_limits<int>::max() to int");
		}
		return static_cast<int>(input);
	}
}

#endif