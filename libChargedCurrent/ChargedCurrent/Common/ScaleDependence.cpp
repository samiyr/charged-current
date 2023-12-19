#ifndef SCALE_DEPENDENCE_H
#define SCALE_DEPENDENCE_H

#include <optional>
#include <concepts>
#include <algorithm>

#include "Common/TRFKinematics.cpp"

// Scale dependence concept. Checks that the object can be called with a TRFKinematics object and that it returns a double.
template <typename T>
concept is_scale_dependence = requires(T scale_function, const TRFKinematics &kinematics) {
	{ scale_function(kinematics) } -> std::same_as<double>;
};

namespace ScaleDependence {
	// A wrapper for a scale dependence function, where the template parameter determines the function type. The function type
	// must satisfy the 'is_scale_dependence' concept; see is_scale_dependence. The function can be left undefined, in which case
	// calling the wrapper defaults to returning the Q^2 value.
	template <is_scale_dependence FunctionType>
	struct Function {
		// Exposes the function type template parameter.
		using type = FunctionType;

		const std::optional<FunctionType> function;
		const double minimum;

		// Constructs an empty function wrapper, in which case calling the constructed object returns Q^2.
		constexpr Function() noexcept : function(std::nullopt), minimum(0.0) {}
		// Constructs a function wrapper, in which case calling the constructed object calls the provided function.
		constexpr Function(const FunctionType func, const double minimum = 0.0) noexcept : function(std::make_optional(func)), minimum(minimum) {}

		// Checks whether the wrapper has a function specified. If false, then calling this object returns Q^2.
		constexpr bool nontrivial() const noexcept {
			return function.has_value();
		}

		// Evaluates the wrapped function if one exists, otherwise returns Q^2.
		constexpr double operator()(const TRFKinematics &kinematics) const noexcept(noexcept(function)) {
			const double value = nontrivial() ? (*function)(kinematics) : kinematics.Q2;
			return std::max(value, minimum);
		}
	};

	// Scale dependence function wrapper that returns Q^2.
	constexpr inline auto trivial = Function([](const TRFKinematics &kinematics) noexcept { return kinematics.Q2; });
	// Scale dependence function wrapper that returns Q^2, with a non-zero minimum.
	[[maybe_unused]] constexpr inline auto clamped(const double minimum) {
		return Function([]([[maybe_unused]] const TRFKinematics &kinematics) noexcept { return kinematics.Q2; }, minimum);
	}
	// Scale dependence function wrapper that returns a constant value.
	[[maybe_unused]] constexpr inline auto constant(const double value) {
		return Function([value]([[maybe_unused]] const TRFKinematics &kinematics) noexcept { return value; });
	}
	// Scale dependence function wrapper that returns Q^2 multiplied by a constant value.
	[[maybe_unused]] constexpr inline auto multiplicative(const double factor) {
		return Function([factor](const TRFKinematics &kinematics) noexcept { return kinematics.Q2 * factor; });
	}
	// Scale dependence function wrapper that returns Q^2 multiplied by a constant value, in such a way that the multiplied value is clamped
	// between a minimum and maximum Q^2 value.
	[[maybe_unused]] constexpr inline auto clamped_multiplicative(const double factor, const double Q2_min, const double Q2_max) {
		return Function([factor, Q2_min, Q2_max](const TRFKinematics &kinematics) noexcept { return std::clamp(kinematics.Q2 * factor, Q2_min, Q2_max); });
	}

	template <is_scale_dependence FunctionType>
	auto operator*(const double lhs, const Function<FunctionType> &rhs) {
		return Function([lhs, &rhs](const TRFKinematics &kinematics) noexcept {
			return lhs * rhs(kinematics);
		}, rhs.minimum);
	}
	
	template <is_scale_dependence FunctionType>
	auto operator*(const Function<FunctionType> &lhs, const double rhs) {
		return Function([&lhs, rhs](const TRFKinematics &kinematics) noexcept {
			return lhs(kinematics) * rhs;
		}, lhs.minimum);
	}
}

// Different options for automatic scale variation of cross sections.
enum class ScaleVariation {
	// No scale variation.
	None, 
	// All applicable scales are varied independently.
	All, 
	// Only the renormalization and factorization scales are varied independently. If a fragmentation scale exists,
	// it's set to the factorization scale.
	RenormalizationFactorization
};

#endif