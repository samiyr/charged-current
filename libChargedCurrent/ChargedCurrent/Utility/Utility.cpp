#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <limits>
#include <algorithm>
#include <ranges>
#include <filesystem>
#include <sstream>
#include <iomanip>

#include "Utility/Math.cpp"

namespace IO {
	inline std::ostream& endl(std::ostream& os) {
		os.put(os.widen('\n'));
		return os;
	}

	bool create_directory_tree(const std::filesystem::path &path) {
		const auto directory = path.parent_path();
		return std::filesystem::create_directories(directory);
	}

	template <typename Integer> requires std::is_integral_v<Integer>
	std::string leading_zeroes(const Integer i, const unsigned int width) {
		std::stringstream ss;
		ss << std::setw(static_cast<int>(width)) << std::setfill('0') << i;
		return ss.str();
	}
}

namespace Collections {
	template <typename T>
	constexpr T intersection(T &v1, T &v2) {
		T v3;

		std::sort(std::begin(v1), std::end(v1));
		std::sort(std::begin(v2), std::end(v2));

		std::set_intersection(std::begin(v1), std::end(v1), std::begin(v2), std::end(v2), std::back_inserter(v3));

		return v3;
	}

	/// Generates n-tuples from a given list. Implicitly assumes that 0 < n <= input.size().
	template <typename T>
	constexpr std::vector<std::vector<T>> tuples(const std::vector<T> input, const std::size_t n) {
		const std::size_t count = Math::ipow(input.size(), n);
		std::vector<std::vector<T>> output(count, std::vector<T>(n));

		for (std::size_t outer = 0; outer < count; outer++) {
			std::size_t current = outer;
			for (std::size_t inner = 0; inner < n; inner++) {
				const std::size_t index = current % input.size();
				output[outer][n - inner - 1] = input[index];
				current /= n;
			}
		}

		return output;
	}

	template <typename T, typename S>
	bool contains(const T &collection, const S &value) {
		return std::find(std::begin(collection), std::end(collection), value) != std::end(collection);
	}

	template <typename T>
	auto closest(const std::vector<T> &input, const T &value) {
		using difference_type = typename std::vector<T>::difference_type;

		const auto iterator = std::lower_bound(std::begin(input), std::end(input), value);
		if (iterator == std::end(input)) { return difference_type{-1}; }

		return std::distance(input.begin(), iterator);
	}
}

template <typename T, typename U>
constexpr T& operator*=(T &v1, U &&v2) {
	for (const auto &[e1, e2] : std::views::zip(v1, std::forward<U>(v2))) {
		e1 *= e2;
	}
	return v1;
}

template<class Specialization, template<typename...> class TemplateClass, typename ...PartialSpecialisation>
concept is_instance = requires (Specialization s) {
    []<typename ...TemplateArgs>(TemplateClass<PartialSpecialisation..., TemplateArgs...>&){}(s);
};

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
	/// Casts a std::size_t to an int with no bounds checking. Supplying a value larger than
	/// std::numeric_limits<int>::max() results in implementation-defined behavior.
	constexpr int size_to_int(const std::size_t input) noexcept {
		return static_cast<int>(input);
	}

	constexpr std::string bool_to_string(const bool input) {
		return input ? "true" : "false";
	}
}

namespace Comparison {
	bool relative_comparison(double x, double y, double epsilon) {
 	   return std::abs(x - y) <= epsilon * std::max(std::abs(x), std::abs(y));
	}
	// bool double_comparison(double a, double b, double tolerance = 1e-5) {
	// 	bool flag = abs(a - b) < tolerance;
	// 	// const bool flag = std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;

	// 	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	// 	if (flag) {
	// 		std::cout << "PASS" << IO::endl;
	// 	} else {
	// 		std::cout << "FAIL" << IO::endl;
	// 	}
	// 	return flag;
	// }
	// bool double_comparison_rel(double a, double b, double tolerance = 1e-5) {
	// 	const bool flag = std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;

	// 	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	// 	if (flag) {
	// 		std::cout << "PASS" << IO::endl;
	// 	} else {
	// 		std::cout << "FAIL" << IO::endl;
	// 	}
	// 	return flag;
	// }

	// double order_of_magnitude(const double x) {
	// 	return std::log10(std::abs(x));
	// }

	// bool double_comparison_digits(double a, double b, int digits = 5) {
	// 	const double difference = order_of_magnitude(a - b);
	// 	const double original = order_of_magnitude(a);
	// 	bool flag = abs(difference - original) >= double(digits);

	// 	std::cout << "Floating-point comparison between " << a << " and " << b << ": ";
	// 	if (flag) {
	// 		std::cout << "PASS" << IO::endl;
	// 	} else {
	// 		std::cout << "FAIL" << IO::endl;
	// 	}
	// 	return flag;
	// }
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

	unsigned int get_default_thread_count() {
		const unsigned int hardware = std::thread::hardware_concurrency();
		return std::max(1U, hardware);
	}
}

#endif