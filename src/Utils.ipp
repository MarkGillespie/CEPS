namespace CEPS {
namespace ImplementationDetails {

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data) {
    out << "{";
    for (size_t i = 0; i < fmin(data.size(), 5); ++i) {
        out << data[i];
        if (i + 1 < data.size()) out << ", ";
    }
    if (data.size() >= 5) out << "... (length " << data.size() << ")";
    out << "}";
    return out;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& data) {
    out << "{";
    for (size_t i = 0; i < fmin(N, 5); ++i) {
        out << data[i];
        if (i + 1 < N) out << ", ";
    }
    if (N >= 5) out << "... (length " << N << ")";
    out << "}";
    return out;
}

template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& data) {
    out << "(" << std::get<0>(data) << ", " << std::get<1>(data) << ")";
    return out;
}

} // namespace ImplementationDetails
} // namespace CEPS
