#include "FlipFormulas.h"
/*
 *          k
 *        / |\
 *       /    \
 *     eki    ejk
 *     /        \
 *   |/          \
 *  i --- eij ----> j
 *    \          /|
 *     \        /
 *    eil      elj
 *       \    /
 *        \| /
 *         l
 */

namespace CEPS {
// Delaunay flip quantity (eq 9) for edge e, and flip(e)
std::pair<double, double>
delaunayCondition(const std::array<double, 5>& lengths) {
    // unpack lengths into conveniently named variables
    double lij, ljk, lki, lil, llj;
    std::tie(lij, ljk, lki, lil, llj) = std::tuple_cat(lengths);

    double llk = flipHyperbolicLength(lengths);

    double flipQ = (lil * lki + ljk * llj) * (lil * ljk + lki * llj) -
                   lij * lij * (ljk * lki + lil * llj);

    double flippedFlipQ = (ljk * lki + lil * llj) * (lil * ljk + lki * llj) -
                          llk * llk * (lil * lki + ljk * llj);
    return std::make_pair(flipQ, flippedFlipQ);
}

size_t flipNormalCoordinate(const std::array<size_t, 5>& normalCoordinates) {
    // unpack normal coordinates into conveniently named variables
    size_t nij, njk, nki, nil, nlj;
    std::tie(nij, njk, nki, nil, nlj) = std::tuple_cat(normalCoordinates);

    int Eilj = positivePart(nlj - nij - nil);
    int Ejil = positivePart(nil - nlj - nij);
    int Elji = positivePart(nij - nil - nlj);

    int Eijk = positivePart(njk - nki - nij);
    int Ejki = positivePart(nki - nij - njk);
    int Ekij = positivePart(nij - njk - nki);

    double Cilj = -(negativePart(nlj - nij - nil) + Ejil + Elji) / 2.;
    double Cjil = -(negativePart(nil - nlj - nij) + Eilj + Elji) / 2.;
    double Clji = -(negativePart(nij - nil - nlj) + Eilj + Ejil) / 2.;

    double Cijk = -(negativePart(njk - nki - nij) + Ejki + Ekij) / 2.;
    double Cjki = -(negativePart(nki - nij - njk) + Eijk + Ekij) / 2.;
    double Ckij = -(negativePart(nij - njk - nki) + Eijk + Ejki) / 2.;

    int doubledAnswer = 2 * Clji + 2 * Ckij + abs(Cjil - Cjki) +
                        abs(Cilj - Cijk) - Elji - Ekij + 2 * Eilj + 2 * Eijk +
                        2 * Ejil + 2 * Ejki - 2 * negativePart(nij - 1);

    return positivePart(std::round(doubledAnswer / 2));
}

double flipEuclideanLength(const std::array<double, 5>& lengths) {
    // unpack lengths into conveniently named variables
    double lij, ljk, lki, lil, llj;
    std::tie(lij, ljk, lki, lil, llj) = std::tuple_cat(lengths);

    // Lay out triangle pair in the plane, and compute the length from k to l
    Vector2 pl{0., 0.};
    Vector2 pj{llj, 0.};
    Vector2 pi = layoutTriangleVertex(
        pl, pj, lij, lil); // involves more arithmetic than strictly necessary
    Vector2 pk = layoutTriangleVertex(pi, pj, ljk, lki);
    double llk = (pl - pk).norm();

    return llk;
}

double flipHyperbolicLength(const std::array<double, 5>& lengths) {
    // unpack lengths into conveniently named variables
    double lij, ljk, lki, lil, llj;
    std::tie(lij, ljk, lki, lil, llj) = std::tuple_cat(lengths);

    return (ljk * lil + lki * llj) / lij;
}

std::array<size_t, 2>
flipRoundabout(const std::array<size_t, 5>& normalCoordinates, size_t nlk,
               size_t dk, size_t dl, size_t rki, size_t rlj) {
    // unpack normal coordinates into conveniently named variables
    size_t nij, njk, nki, nil, nlj;
    std::tie(nij, njk, nki, nil, nlj) = std::tuple_cat(normalCoordinates);

    /*         Rotates edge clockwise
     *          k                   k
     *        / |\                / ||\
     *       /    \              /  |  \
     *     eki     ejk          /   |   \
     *     /        \         e2    |    e1
     *   |/          \       |/     |     \
     *  i --- eij ---> j    i \     v     /| j
     *    \          /|        \    |    /
     *     \        /           \   |   /
     *    eil      elj          e3  |  e4
     *       \    /               \ | /
     *        \| /                _\|/
     *         l                    l
     */

    size_t deltaki = (nki == 0) ? 1 : 0;
    size_t deltalj = (nlj == 0) ? 1 : 0;

    size_t ekil = positivePart(nil - nki - nlk);
    size_t eljk = positivePart(njk - nlj - nlk);

    size_t rkl = (rki + deltaki + ekil) % dk;
    size_t rlk = (rlj + deltalj + eljk) % dl;

    return {rkl, rlk};
}
} // namespace CEPS
