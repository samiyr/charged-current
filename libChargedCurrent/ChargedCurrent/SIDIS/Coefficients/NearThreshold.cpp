#ifndef SIDIS_FUNCTIONS_NEAR_THRESHOLD_H
#define SIDIS_FUNCTIONS_NEAR_THRESHOLD_H

#include <cstddef>
#include <cmath>

namespace SIDISFunctions::NearThreshold {
	constexpr double d(const double x, const double xq) {
		return xq / (1.0 - x);
	}
	constexpr double D(const std::size_t n, const double x, const double xi, const double xq, const double xq_hat) {
		const double difference = xq_hat - xq;
		const double log_term_1 = std::pow(std::log(1.0 - x), 1 + n) / (1.0 + static_cast<double>(n));
		const double log_term_2 = std::pow(std::log(1.0 - xi), n) / (1.0 - xi);

		return difference * log_term_2 + xq * log_term_1 / (1.0 - x);
	}
	constexpr double l(const std::size_t n, const double xi, const double xq_hat) {
		return xq_hat * std::pow(std::log(1.0 - xi), n);
	}

}

#endif