#ifndef PDF_CONCEPT_H
#define PDF_CONCEPT_H

#include <concepts>
#include <string>

#include "Common/Flavor.cpp"

#include "PDF/Interfaces/LHAInterface.cpp"
#include "PDF/Interfaces/LHASetInterface.cpp"
#include "PDF/Interfaces/FunctionalFormInterface.cpp"

template <typename T>
concept is_pdf_interface = requires(T pdf, const FlavorType flavor, const double x, const double Q2) {
    { pdf.evaluate(x, Q2) } -> std::same_as<void>;
    { pdf.xf_evaluate(flavor, x, Q2) } -> std::same_as<double>;
    { pdf.xf(flavor) } -> std::same_as<double>;
    { pdf.xg() } -> std::same_as<double>;
    { pdf.alpha_s(Q2) } -> std::same_as<double>;
	{ pdf.quark_mass(flavor) } -> std::same_as<double>;
};

static_assert(is_pdf_interface<LHAInterface<>>, "LHAInterface doesn't satisfy is_pdf_interface");
static_assert(
	is_pdf_interface<FunctionalFormInterface<decltype([](const FlavorType, const double, const double) { return 0.0; })>>, 
	"FunctionalFormInterface doesn't satisfy is_pdf_interface"
);

template <typename T>
concept is_pdf_set_interface = requires(const T pdf, const FlavorType flavor, const typename T::size_type index) {
	{ pdf[index] };
	requires is_pdf_interface<decltype(pdf[index])>;

	{ pdf.central() };
	requires is_pdf_interface<decltype(pdf.central())>;

	{ pdf.quark_mass(flavor) } -> std::same_as<double>;
	{ pdf.size() } -> std::same_as<typename T::size_type>;
};

static_assert(is_pdf_set_interface<LHASetInterface<>>, "LHASetInterface doesn't satisfy is_pdf_set_interface");

template<typename T>
concept has_pdf_activation = requires(const T &pdf) {
    { pdf.activate() } -> std::same_as<void>;
};

#endif