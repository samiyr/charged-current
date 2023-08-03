#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <limits>
#include <algorithm>
#include <ranges>

template <typename T>
static constexpr T POW2(const T x) {
	return x * x;
}

template <typename T>
static constexpr T POW4(const T x) {
	return x * x * x * x;
}

namespace IO {
	inline std::ostream& endl(std::ostream& os) {
		os.put(os.widen('\n'));
		return os;
	}
}

namespace Collections {
	template <typename T>
	constexpr std::vector<T> vector_intersection(std::vector<T> &v1, std::vector<T> &v2) {
		std::vector<T> v3;

		std::sort(v1.begin(), v1.end());
		std::sort(v2.begin(), v2.end());

		std::set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), std::back_inserter(v3));

		return v3;
	}

	template <typename T>
	static void multiply(std::vector<T> &v1, const std::vector<T> &v2) {
		for (typename std::vector<T>::size_type i = 0; i < v1.size(); i++) {
			v1[i] *= v2[i];
		}
	}
}

template <typename Type, std::size_t Size>
constexpr std::array<Type, Size> operator*(const std::array<Type, Size> &lhs, const std::array<Type, Size> &rhs) {
	std::array<Type, Size> output{};

	for (const auto &[lhs_value, rhs_value, output_value] : std::views::zip(lhs, rhs, output)) {
		output_value = lhs_value * rhs_value;
	}

	return output;
}

template <typename T>
constexpr std::vector<T> operator+(const std::vector<T> &lhs, const T &rhs) {
	std::vector<T> output(lhs.size());
	std::transform(lhs.begin(), lhs.end(), output.begin(), [&rhs](const T &current) { return current + rhs; });
	return output;
}

template <typename T>
constexpr std::vector<T> operator-(const std::vector<T> &lhs, const T &rhs) {
	std::vector<T> output(lhs.size());
	std::transform(lhs.begin(), lhs.end(), output.begin(), [&rhs](const T &current) { return current - rhs; });
	return output;
}

namespace Conversion {
	constexpr int size_to_int(const size_t input) {
		if (input > std::numeric_limits<int>::max()) {
			throw std::runtime_error("Cannot cast size_t larger than std::numeric_limits<int>::max() to int");
		}
		return static_cast<int>(input);
	}

	constexpr std::string bool_to_string(const bool input) {
		return input ? "true" : "false";
	}
}

namespace Comparison {
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
}

namespace Utility {
	void non_aborting_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
		std::cout << "gsl error " << gsl_errno << " (" << std::string(file) << ", " << line << "): " << std::string(reason) << IO::endl;
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
}

#endif