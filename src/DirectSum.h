#pragma once

#include <sstream>

namespace CEPS {
namespace ImplementationDetails {

template <typename T1, typename T2>
struct DirectSum {
    T1 x;
    T2 y;
    DirectSum(){};
    DirectSum(T1 x_, T2 y_) : x(x_), y(y_){};

    double norm() const;
    double norm2() const;
};

template <typename T1, typename T2>
DirectSum<T1, T2> operator+(const DirectSum<T1, T2>& a,
                            const DirectSum<T1, T2>& b);

template <typename T1, typename T2>
DirectSum<T1, T2> operator-(const DirectSum<T1, T2>& a,
                            const DirectSum<T1, T2>& b);

template <typename T1, typename T2>
DirectSum<T1, T2> operator*(const double& s, const DirectSum<T1, T2>& a);

template <typename T1, typename T2>
DirectSum<T1, T2> operator*(const DirectSum<T1, T2>& a, const double& s);

template <typename T1, typename T2>
DirectSum<T1, T2> operator*(const DirectSum<T1, T2>& a, const double& s);

double dot(double a, double b);

template <typename T1, typename T2>
double dot(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b);

template <typename T1, typename T2>
bool operator==(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b);

template <typename T1, typename T2>
bool operator!=(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b);

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& str, const DirectSum<T1, T2>& a);

} // namespace ImplementationDetails
} // namespace CEPS

#include "DirectSum.ipp"
