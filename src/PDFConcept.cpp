#ifndef PDF_CONCEPT_H
#define PDF_CONCEPT_H

#include <concepts>
#include <string>
#include "Flavor.cpp"

template <typename T>
concept PDFConcept = requires(T pdf, const FlavorType flavor, const double x, const double Q2) {
    { pdf.evaluate(x, Q2) } -> std::same_as<void>;
    { pdf.xf_evaluate(flavor, x, Q2) } -> std::same_as<double>;
    { pdf.xf(flavor) } -> std::same_as<double>;
    { pdf.xg() } -> std::same_as<double>;
    { pdf.alpha_s(Q2) } -> std::same_as<double>;
};

#endif