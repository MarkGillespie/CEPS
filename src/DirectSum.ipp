namespace CEPS {
namespace ImplementationDetails {
template <typename T1, typename T2>
double DirectSum<T1, T2>::norm() const {
    return sqrt(norm2());
}

template <typename T1, typename T2>
double DirectSum<T1, T2>::norm2() const {
    return dot(*this, *this);
}

template <typename T1, typename T2>
DirectSum<T1, T2> operator+(const DirectSum<T1, T2>& a,
                            const DirectSum<T1, T2>& b) {
    return DirectSum<T1, T2>(a.x + b.x, a.y + b.y);
}

template <typename T1, typename T2>
DirectSum<T1, T2> operator-(const DirectSum<T1, T2>& a,
                            const DirectSum<T1, T2>& b) {
    return DirectSum<T1, T2>(a.x - b.x, a.y - b.y);
}

template <typename T1, typename T2>
DirectSum<T1, T2> operator*(const double& s, const DirectSum<T1, T2>& a) {
    return DirectSum<T1, T2>(s * a.x, s * a.y);
}

template <typename T1, typename T2>
DirectSum<T1, T2> operator*(const DirectSum<T1, T2>& a, const double& s) {
    return DirectSum<T1, T2>(a.x * s, a.y * s);
}

template <typename T1, typename T2>
DirectSum<T1, T2> operator/(const DirectSum<T1, T2>& a, const double& s) {
    return DirectSum<T1, T2>(a.x / s, a.y / s);
}

template <typename T1, typename T2>
double dot(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b) {
    return dot(a.x, b.x) + dot(a.y, b.y);
}

template <typename T1, typename T2>
bool operator==(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b) {
    return a.x == b.x && a.y == b.y;
}

template <typename T1, typename T2>
bool operator!=(const DirectSum<T1, T2>& a, const DirectSum<T1, T2>& b) {
    return !(a == b);
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& str, const DirectSum<T1, T2>& a) {
    str << "{" << a.x << ", " << a.y << "}";
    return str;
}
} // namespace ImplementationDetails
} // namespace CEPS
