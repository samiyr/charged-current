#ifndef SCALE_DEPENDENCE_H
#define SCALE_DEPENDENCE_H

#include "TRFKinematics.cpp"
#include <optional>
#include <concepts>

namespace ScaleDependence {
	template <typename T>
	concept Concept = requires(T scale_function, const TRFKinematics &kinematics) {
		{ scale_function(kinematics) } -> std::same_as<double>;
	};

	template <Concept FunctionType>
	struct Function {
		using type = FunctionType;

		constexpr Function() : function(std::nullopt) {}
		constexpr Function(const FunctionType func) : function(std::make_optional(func)) {}

		constexpr bool nontrivial() const {
			return function.has_value();
		}

		constexpr double compute(const TRFKinematics &kinematics) const {
			if (nontrivial()) {
				return (*function)(kinematics);
			}
			return kinematics.Q2;
		}

		constexpr double operator()(const TRFKinematics &kinematics) const {
			return compute(kinematics);
		}

		private:
		const std::optional<FunctionType> function;
	};

	constexpr static auto trivial = Function([](const TRFKinematics &kinematics) { return kinematics.Q2; });
	[[maybe_unused]] constexpr static auto constant(const double value) {
		return Function([value]([[maybe_unused]] const TRFKinematics &kinematics) { return value; });
	}
	[[maybe_unused]] constexpr static auto multiplicative(const double factor) {
		return Function([factor](const TRFKinematics &kinematics) { return kinematics.Q2 * factor; });
	}
	// [[maybe_unused]] static auto Q2 = std::make_optional([](const TRFKinematics &kinematics) { return kinematics.Q2; });
	// static auto trivial = [](const TRFKinematics &kinematics) { return kinematics.Q2; };
	// [[maybe_unused]] static auto constant(const double value) {
	// 	return std::make_optional([value]([[maybe_unused]] const TRFKinematics &kinematics) { return value; });
	// }
	// [[maybe_unused]] static auto multiplicative(const double factor) {
	// 	return std::make_optional([factor](const TRFKinematics &kinematics) { return kinematics.Q2 * factor; });
	// }
}

#endif