#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include <cmath>

namespace Coefficients {
	namespace DIS {
		double C2q(double x) {
			double a = (1 + x * x) / (1 - x);
			double b = std::log1p(-x) - std::log(x);
			double c = 0.25 * (9 + 5 * x);

			double value = a * (b - 3.0 / 4.0) + c;

			return (4.0 / 3.0) * value;
		}

		double C2g(double x) {
			double a = x * x + (1 - x) * (1 - x);
			double b = std::log1p(-x) - std::log(x);
			double c = 4 * x * (1 - x);

			return 0.5 * a * b - 0.5 + c;
		}
	}
}

#endif