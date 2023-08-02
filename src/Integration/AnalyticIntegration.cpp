// #ifndef ANALYTIC_INTEGRATION_H
// #define ANALYTIC_INTEGRATION_H

// #include <string>

// template <auto Value>
// struct NumericalConstant {
// 	static constexpr auto get = Value;
// };

// // https://dev.to/sgf4/strings-as-template-parameters-c20-4joh
// template<std::size_t N>
// struct StringLiteral {
//     char data[N] {};

//     consteval StringLiteral(const char (&str)[N]) {
//         std::copy_n(str, N, data);
//     }

//     consteval bool operator==(const StringLiteral<N> str) const {
//         return std::equal(str.data, str.data+N, data);
//     }

//     template<std::size_t N2>
//     consteval bool operator==(const StringLiteral<N2> s) const {
//         return false;
//     }

//     template<std::size_t N2>
//     consteval StringLiteral<N+N2-1> operator+(const StringLiteral<N2> str) const {
//         char newchar[N+N2-1] {};
//         std::copy_n(data, N-1, newchar);
//         std::copy_n(str.data, N2, newchar+N-1);
//         return newchar;
//     }

//     consteval char operator[](std::size_t n) const {
//         return data[n];
//     }

//     consteval std::size_t size() const {
//         return N-1;
//     }
// };


// template <StringLiteral Name>
// struct Variable {
// 	static constexpr auto evaluate = [](const auto x) {
// 		return x;
// 	}
// };

// template <auto Coefficient, auto Exponent, StringLiteral Variable>
// struct Monomial {
// 	static constexpr auto evaluate = [](const auto x) {
// 		return Coefficient * std::pow(x, Exponent);
// 	}
// };

// template <typename Left, typename Right>
// struct Plus {
// 	// static constexpr auto evaluate = [](const auto x) {
// 	// 	return Left::evaluate(x) + Right::evaluate(x);
// 	// }
// };

// template <typename Left, typename Right>
// struct Times;

// template <auto Center>
// struct DiracDelta;

// template <typename Integrand, StringLiteral Variable>
// struct Integral;

// template <auto Coefficient, auto Exponent, StringLiteral Variable>
// struct Integral<Monomial<Coefficient, Exponent, Variable>, Variable> {	
// 	static constexpr auto integrate = [](const auto lower, const auto upper) {
// 		if constexpr (Exponent == -1) {
// 			return Coefficient * std::log(upper / lower);
// 		}
// 		return Coefficient * (std::pow(upper, 1 + Exponent) - std::pow(lower, 1 + Exponent)) / (1 + Exponent);
// 	};

// };

// template <StringLiteral Name>
// struct Integral<Variable<Name>, Name> {
// 	static constexpr auto integrate = Integral<Monomial<1, 1, Name>, Name>::integrate;
// };

// template <typename Left, typename Right, StringLiteral Variable>
// struct Integral<Plus<Left, Right>, Variable> {
// 	static constexpr auto evaluator = [](const auto lower, const auto upper) {
// 		return Integral<Left, Variable>::evaluator(lower, upper) + Integral<Right, Variable>::evaluator(lower, upper);
// 	};
// };

// template <auto Center, typename Right, StringLiteral Variable>
// struct Integral<Times<DiracDelta<Center>, Right>> {
// 	static constexpr auto evaluator = [](const auto lower, const auto upper) {

// 	}
// }

// #endif