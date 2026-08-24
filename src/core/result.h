#ifndef MISPITOOLS_CORE_RESULT_H
#define MISPITOOLS_CORE_RESULT_H

#include <optional>
#include <string>
#include <utility>

namespace mispitools {
namespace core {

/// @brief Carrier for a possibly-failing computation in the kernel.
///
/// `value` is populated iff `error` is empty. The kernel never throws
/// across its boundary (see DESIGN.md §5); errors flow back as values
/// and are translated by the binding layer (Rcpp::stop, JS exceptions).
template <typename T>
struct Result {
    std::optional<T> value;
    std::string error;

    bool ok() const noexcept { return value.has_value(); }
    const T& operator*() const { return *value; }
    T& operator*() { return *value; }
};

template <typename T>
inline Result<T> ok_result(T v) {
    return Result<T>{ std::optional<T>{ std::move(v) }, std::string{} };
}

template <typename T>
inline Result<T> err_result(std::string msg) {
    return Result<T>{ std::nullopt, std::move(msg) };
}

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_RESULT_H
