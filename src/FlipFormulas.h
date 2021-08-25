#pragma once

#include "Utils.h"

#include "geometrycentral/utilities/elementary_geometry.h"

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
// Delaunay flip quantity (eq ???) for edge e, and flip(e)
std::pair<double, double>
delaunayCondition(const std::array<double, 5>& lengths);

size_t flipNormalCoordinate(const std::array<size_t, 5>& normalCoordinates);

double flipEuclideanLength(const std::array<double, 5>& lengths);

double flipHyperbolicLength(const std::array<double, 5>& lengths);

std::array<size_t, 2>
flipRoundabout(const std::array<size_t, 5>& normalCoordinates, size_t nlk,
               size_t dk, size_t dl, size_t rki, size_t rlj);
} // namespace CEPS
