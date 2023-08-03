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

		constexpr Function() noexcept : function(std::nullopt) {}
		constexpr Function(const FunctionType func) noexcept : function(std::make_optional(func)) {}

		constexpr bool nontrivial() const noexcept {
			return function.has_value();
		}

		constexpr double operator()(const TRFKinematics &kinematics) const noexcept(noexcept(function)) {
			if (nontrivial()) {
				return (*function)(kinematics);
			}
			return kinematics.Q2;		
		}

		private:
		const std::optional<FunctionType> function;
	};

	constexpr static auto trivial = Function([](const TRFKinematics &kinematics) noexcept { return kinematics.Q2; });
	[[maybe_unused]] constexpr static auto constant(const double value) {
		return Function([value]([[maybe_unused]] const TRFKinematics &kinematics) noexcept { return value; });
	}
	[[maybe_unused]] constexpr static auto multiplicative(const double factor) {
		return Function([factor](const TRFKinematics &kinematics) noexcept { return kinematics.Q2 * factor; });
	}
}

#endif