#ifndef MISPITOOLS_CORE_MARKER_H
#define MISPITOOLS_CORE_MARKER_H

#include <cstdint>
#include <string>
#include <vector>

namespace mispitools {
namespace core {

using AlleleIndex = std::int32_t;
using GenotypeIndex = std::int32_t;
constexpr AlleleIndex kNoAllele = -1;

/// @brief Single forensic marker with population allele frequencies.
///
/// `numeric_labels` parallels `freqs` and carries the allele labels as
/// doubles for stepwise mutation (NaN entries are accepted but trigger
/// an error if the stepwise builder later consumes them).
struct Marker {
    std::string id;
    AlleleIndex n_alleles = 0;
    std::vector<double> freqs;            // size n_alleles, sums to 1 within tol
    std::vector<double> numeric_labels;   // size n_alleles, used by stepwise
};

/// @brief Number of unordered genotypes G = K*(K+1)/2.
inline GenotypeIndex n_genotypes(AlleleIndex K) noexcept {
    return static_cast<GenotypeIndex>(K) * (K + 1) / 2;
}

/// @brief Map an unordered allele pair (0 <= a1 <= a2 < K) to its genotype
/// index in the column-major-lex order used by the R reference engine.
///   g(a1, a2) = a2 * (a2 + 1) / 2 + a1.
/// For K=3 the enumeration is (0,0), (0,1), (1,1), (0,2), (1,2), (2,2).
inline GenotypeIndex pair_to_genotype_index(
        AlleleIndex a1, AlleleIndex a2) noexcept {
    return static_cast<GenotypeIndex>(a2) * (a2 + 1) / 2 + a1;
}

/// @brief Inverse of pair_to_genotype_index().
inline void genotype_index_to_pair(
        GenotypeIndex g, AlleleIndex& a1, AlleleIndex& a2) noexcept {
    AlleleIndex a = 0;
    while (static_cast<GenotypeIndex>(a + 1) * (a + 2) / 2 <= g) {
        ++a;
    }
    a2 = a;
    a1 = static_cast<AlleleIndex>(
        g - static_cast<GenotypeIndex>(a2) * (a2 + 1) / 2);
}

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int marker_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_MARKER_H
