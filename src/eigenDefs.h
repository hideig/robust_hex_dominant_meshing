#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <queue>
#include <set>
#include <tuple>
#include <vector>

using IndexType = long long;

template<typename T>
using SMatrix = Eigen::SparseMatrix<T, Eigen::RowMajor>;

template<typename T>
using Triple = Eigen::Triplet<T>;
using Tripled = Triple<double>;

using _Matrix_3 = double[3][3];

using Matrix_3 = Eigen::Matrix<double, 3, 3>;

using Triangle = Eigen::Matrix<IndexType, 3, 1>;

template<int n>
using FaceVertices = std::array<IndexType, n>;

using V3d = Eigen::Vector3d;

using V4d = Eigen::Vector4d;

template<typename T>
using ColMatrix =
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template<int columns>
using RowVectors =
  Eigen::Matrix<double, Eigen::Dynamic, columns, Eigen::RowMajor>;

constexpr double PI = 3.141592653589793238463;

using IndexPair = std::array<IndexType, 2>;