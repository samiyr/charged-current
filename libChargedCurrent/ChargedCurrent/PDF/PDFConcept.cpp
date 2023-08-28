#ifndef PDF_CONCEPT_H
#define PDF_CONCEPT_H

#include <concepts>
#include <string>

#include "Common/Flavor.cpp"

template <typename T>
concept is_pdf_interface = requires(T pdf, const FlavorType flavor, const double x, const double Q2) {
    { pdf.evaluate(x, Q2) } -> std::same_as<void>;
    { pdf.xf_evaluate(flavor, x, Q2) } -> std::same_as<double>;
    { pdf.xf(flavor) } -> std::same_as<double>;
    { pdf.xg() } -> std::same_as<double>;
    { pdf.alpha_s(Q2) } -> std::same_as<double>;
};

template<typename T>
concept has_pdf_activation = requires(const T &pdf) {
    { pdf.activate() } -> std::same_as<void>;
};

#endif