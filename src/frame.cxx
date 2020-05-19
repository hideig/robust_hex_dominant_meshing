#include "frame.h"
#include "iostream"

namespace LatticeCore {
void
S6::constructor(int a, int b, int c)
{
  trans_matrix_ = Eigen::Matrix<int, 3, 3>::Zero();
  trans_matrix_(0, a) = 1;
  trans_matrix_(1, b) = 1;
  trans_matrix_(2, c) = 1;
}

const S6 S6::TRANS_MATRICES[6] = { S6(0), S6(1), S6(2), S6(3), S6(4), S6(5) };

S6::S6(int index)
{
  constructor(Permutation_3::permutations[index][0],
              Permutation_3::permutations[index][1],
              Permutation_3::permutations[index][2]);
}
S6::S6(int a, int b, int c)
{
  constructor(a, b, c);
}

Eigen::Matrix<int, 3, 1> S6::operator*(Eigen::Matrix<int, 3, 1>& rv)
{
  return trans_matrix_ * rv;
}
int
S6::index_transform(int u, int v)
{
  auto U = TRANS_MATRICES[u];
  auto V = TRANS_MATRICES[v];
  auto Z = V.trans_matrix_ * (U.trans_matrix_);

  std::array<int, 3> p;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (Z(i, j) == 1) {
        p[j] = i;
        break;
      }
    }
  }
  for (int i = 0; i < 6; ++i) {
    if (p == Permutation_3::permutations[i].data_) {
      return i;
    }
  }
  return -1;
}

const Permutation_3 Permutation_3::permutations[6] = {
  Permutation_3(0, 1, 2), Permutation_3(0, 2, 1), Permutation_3(1, 0, 2),
  Permutation_3(1, 2, 0), Permutation_3(2, 0, 1), Permutation_3(2, 1, 0)
};

Permutation_3 Permutation_3::operator*(const Permutation_3& rv) const
{
  int res[3];
  for (int i = 0; i < 3; ++i) {
    res[i] = rv.data_[data_[i]];
  }
  return Permutation_3(res[0], res[1], res[2]);
}
bool
Permutation_3::operator==(const Permutation_3& rv) const
{
  return data_ == rv.data_;
}
int
Permutation_3::transform(int lv, int rv)
{
  Permutation_3 lp = permutations[lv];
  Permutation_3 rp = permutations[rv];
  Permutation_3 ans = lp * rp;
  for (int i = 0; i < 6; ++i) {
    if (permutations[i] == ans) {
      return i;
    }
  }
  return -1;
}

double
Frame3d::harmonic(double x, double y, double z)
{
  return x * x * y * y + y * y * z * z + z * z * x * x;
}

double
Frame3d::harmonic(const V3d& x)
{
  return x(1) * x(1) * x(2) * x(2) + x(0) * x(0) * x(2) * x(2) +
         x(0) * x(0) * x(2) * x(2);
}

double
Frame3d::H(const Matrix_3& Rst)
{
  double res = 0.0;
  for (size_t i = 0; i < 3; i++) {
    res += harmonic(Rst(i, 0), Rst(i, 1), Rst(i, 2));
  }
  return res;
}

Matrix_3
Frame3d::der_H(const Matrix_3& Rst)
{
  Matrix_3 res;
  for (size_t i = 0; i < 3; i++) {
    res(i, 0) = harmonic_x(Rst(i, 0), Rst(i, 1), Rst(i, 2));
    res(i, 1) = harmonic_y(Rst(i, 0), Rst(i, 1), Rst(i, 2));
    res(i, 2) = harmonic_z(Rst(i, 0), Rst(i, 1), Rst(i, 2));
  }
  return res;
}

double
Frame3d::harmonic_x(double x, double y, double z)
{
  return x * y * y * 2 + x * z * z * 2;
}
double
Frame3d::harmonic_x(const V3d& u)
{
  return u(0) * u(1) * u(1) * 2 + u(0) * u(2) * u(2) * 2;
}

double
Frame3d::harmonic_y(double x, double y, double z)
{
  return y * x * x * 2 + y * z * z * 2;
}
double
Frame3d::harmonic_y(const V3d& u)
{
  return u(1) * u(0) * u(0) * 2 + u(1) * u(3) * u(3) * 2;
}

double
Frame3d::harmonic_z(double x, double y, double z)
{
  return z * x * x * 2 + z * y * y * 2;
}
double
Frame3d::harmonic_z(const V3d& u)
{
  return u(2) * u(0) * u(0) * 2 + u(2) * u(1) * u(1) * 2;
}

double
Frame3d::metric(const Matrix_3& A, const Matrix_3& B)
{
  Matrix_3 Rab = A.transpose() * B;
  Matrix_3 Rba = B.transpose() * A; // Rab^T = Rba
  return H(Rab) + H(Rba);
}

Frame3d
Frame3d::average(const Frame3d& a, const Frame3d& b)
{

  return Frame3d();
}

double
Frame3d::diff(const Frame3d& b, const Permutation_3& permute) const
{

  double dff = 0;
  /* compare major principal vectors */

  dff += sin(m_.col(0), b.m_.col(permute[0]));

  /* compare mid pricipal vectors */

  dff += sin(m_.col(1), b.m_.col(permute[1]));

  /* compare minor pricipal vectors */

  dff += sin(m_.col(2), b.m_.col(permute[2]));

  return dff;
}

int
Frame3d::get_S6_matching_index(const Frame3d& b) const
{
  double min_diff = 1e20;
  double tmp;
  int index = -1;
  for (int i = 0; i < 6; ++i) {
    tmp = this->diff(b, Permutation_3::permutations[i]);
    if (tmp < min_diff) {
      min_diff = tmp;
      index = i;
    }
  }
  return index;
}

void
Frame3d::normalize()
{
  for (size_t i = 0; i < 3; i++) {
    m_.col(i) /= m_.col(i).norm();
  }
}

std::ostream&
operator<<(std::ostream& out, const Frame3d& frame)
{
  out << frame.m_(0, 0) << ',' << frame.m_(1, 0) << ',' << frame.m_(2, 0)
      << ',';
  out << frame.m_(0, 1) << ',' << frame.m_(1, 1) << ',' << frame.m_(2, 1)
      << ',';
  out << frame.m_(0, 2) << ',' << frame.m_(1, 2) << ',' << frame.m_(2, 2);
  return out;
}

EulerAngleFrame::EulerAngleFrame()
  : u_(0.0)
  , v_(0.0)
  , w_(0.0)
{}

EulerAngleFrame::EulerAngleFrame(double u, double v, double w)
  : u_(u)
  , v_(v)
  , w_(w)
{}

EulerAngleFrame::EulerAngleFrame(const Matrix_3& r)
{
  u_ = std::atan2(r(2, 1), r(2, 2));
  v_ = std::atan2(-r(2, 0), std::sqrt(r(2, 1) * r(2, 1) + r(2, 2) * r(2, 2)));
  w_ = std::atan2(r(1, 0), r(0, 0));
}

Matrix_3
EulerAngleFrame::get_matrix()
{
  return EulerAngleFrame::R(u_, v_, w_);
}

Matrix_3
EulerAngleFrame::R(double u, double v, double w)
{
  Matrix_3 ret;
  double sinu = std::sin(u);
  double cosu = std::cos(u);
  double sinv = std::sin(v);
  double cosv = std::cos(v);
  double sinw = std::sin(w);
  double cosw = std::cos(w);

  ret << cosv * cosw, cosw * sinu * sinv - cosu * sinw,
    sinu * sinw + cosu * cosw * sinv, cosv * sinw,
    cosu * cosw + sinu * sinv * sinw, cosu * sinv * sinw - cosw * sinu, -sinv,
    cosv * sinu, cosu * cosv;

  return ret;
}
Matrix_3
EulerAngleFrame::dR_1(double u, double v, double w)
{
  Matrix_3 ret;
  double sinu = std::sin(u);
  double cosu = std::cos(u);
  double sinv = std::sin(v);
  double cosv = std::cos(v);
  double sinw = std::sin(w);
  double cosw = std::cos(w);

  ret << 0, sinu * sinw + cosu * cosw * sinv, cosu * sinw - cosw * sinu * sinv,
    0, cosu * sinv * sinw - cosw * sinu, -cosu * cosw - sinu * sinv * sinw, 0,
    cosu * cosv, -cosv * sinu;

  return ret;
}
Matrix_3
EulerAngleFrame::dR_2(double u, double v, double w)
{
  Matrix_3 ret;
  double sinu = std::sin(u);
  double cosu = std::cos(u);
  double sinv = std::sin(v);
  double cosv = std::cos(v);
  double sinw = std::sin(w);
  double cosw = std::cos(w);

  ret << -cosw * sinv, cosv * cosw * sinu, cosu * cosv * cosw, -sinv * sinw,
    cosv * sinu * sinw, cosu * cosv * sinw, -cosv, -sinu * sinv, -cosu * sinv;

  return ret;
}
Matrix_3
EulerAngleFrame::dR_3(double u, double v, double w)
{
  Matrix_3 ret;
  double sinu = std::sin(u);
  double cosu = std::cos(u);
  double sinv = std::sin(v);
  double cosv = std::cos(v);
  double sinw = std::sin(w);
  double cosw = std::cos(w);

  ret << -cosv * sinw, -cosu * cosw - sinu * sinv * sinw,
    cosw * sinu - cosu * sinv * sinw, cosv * cosw,
    cosw * sinu * sinv - cosu * sinw, sinu * sinw + cosu * cosw * sinv, 0, 0, 0;

  return ret;
}

QuaternionFrame
QuaternionFrame::complement() const
{
  auto& u = q.coeffs();
  return QuaternionFrame(-u(0), -u(1), -u(2), u(3));
}

Matrix_3
QuaternionFrame::to_rotation_matrix() const
{
  return q.normalized().toRotationMatrix();
  // return q.toRotationMatrix();
}

void
QuaternionFrame::print() const
{
  double* x = coeffs().data();
  for (int i = 0; i < 4; ++i)
    std::cerr << x[i] << ' ';
  std::cerr << std::endl;
}

QuaternionFrame
QuaternionFrame::average(const QuaternionFrame& a, const QuaternionFrame& b)
{
  return QuaternionFrame(a.q.slerp(0.50, b.rotate_to(a).q));
}

inline Eigen::Vector4d
quaternion_average(const std::vector<Eigen::Vector4d>& q,
                   const std::vector<double>& weights)
{
  double total_weight = 0.0;
  int n = q.size();
  for (double w : weights)
    total_weight += w;
  /******** simple average ********************/
  Eigen::Vector4d aq(0, 0, 0, 0);
  for (int i = 0; i < n; ++i) {
    // aq += list[i].rotate_to(ref).coeffs() * weight[i];
    aq += q[i].normalized() * weights[i];
  }
  aq /= total_weight;
  return aq;

  /********** eigen **********************/
  // Eigen::Matrix4d T = Eigen::Matrix4d::Zero();
  // for (int i = 0; i < n; ++i) {
  //   auto& x = q[i].normalized();
  //   T += x * x.transpose() * weights[i];
  // }
  // T /= total_weight;
  // /* eigen vector according to max eigen value is the averaged quaternion */
  // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 4, 4>> solver(T);
  // auto eig_values = solver.eigenvalues();
  // auto eig_vectors = solver.eigenvectors();
  // int id = 0;
  // double max = eig_values[0];
  // for (int i = 1; i < 4; ++i) {
  //   if (eig_values[i] > max) {
  //     max = eig_values[i];
  //     id = i;
  //   }
  // }
  // return eig_vectors.col(id);
}

QuaternionFrame
QuaternionFrame::weighted_average(const std::vector<QuaternionFrame>& list,
                                  const std::vector<double>& weight,
                                  bool rotate_invariant)
{
  size_t n = list.size();

  std::vector<Eigen::Vector4d> coeffs;
  if (rotate_invariant) {
    const auto& ref = list[0];
    coeffs.push_back(ref.coeffs());
    for (int i = 1; i < n; ++i) {
      coeffs.push_back(list[i].rotate_to(ref).coeffs());
    }
  } else {
    for (auto& q : list) {
      coeffs.push_back(q.coeffs());
    }
  }
  Eigen::Vector4d aq = quaternion_average(coeffs, weight);
  QuaternionFrame q(aq(0), aq(1), aq(2), aq(3));
  q.normalize();
  return q;
}

/**
 * @brief find a transformation that make q1,q2 as close as possible
 *
 * @param q1
 * @param q2
 * @return QuaternionFrame r, which minimize d(q1,q2*r)
 */
QuaternionFrame
QuaternionFrame::request_best_match(const QuaternionFrame& q1,
                                    const QuaternionFrame& q2)
{
  // use Gao Xifeng's code
  double dp[4] = { q1.dot(QuaternionFrame(q2.w(), q2.z(), -q2.y(), -q2.x())),
                   q1.dot(QuaternionFrame(-q2.z(), q2.w(), q2.x(), -q2.y())),
                   q1.dot(QuaternionFrame(q2.y(), -q2.x(), q2.w(), -q2.z())),
                   q1.dot(QuaternionFrame(q2.x(), q2.y(), q2.z(), q2.w())) };

  double a[4] = {
    std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]), std::abs(dp[3])
  };

  // find M, N, max indices of a[i]
  int M = 0, N = 1, m = 2, n = 3;
  if (a[M] < a[m])
    std::swap(M, m);
  if (a[N] < a[n])
    std::swap(N, n);
  if (a[M] < a[N]) {
    std::swap(M, N);
    m = n;
  }
  if (a[N] < a[m])
    N = m;

  const double s = std::sqrt(.5f);
  const double h = .5f;

  double vA = a[M];
  double vB = (a[M] + a[N]) * s;
  double vC = (a[0] + a[1] + a[2] + a[3]) * h;

  QuaternionFrame result(0, 0, 0, 0);
  if (vA > vB && vA > vC) {
    result[M] = std::copysign(1.f, dp[M]);
  } else if (vB > vC) {
    result[M] = std::copysign(s, dp[M]);
    result[N] = std::copysign(s, dp[N]);
  } else {
    result[0] = std::copysign(h, dp[0]);
    result[1] = std::copysign(h, dp[1]);
    result[2] = std::copysign(h, dp[2]);
    result[3] = std::copysign(h, dp[3]);
  }

  return result;
}

QuaternionFrame
QuaternionFrame::rotate_to(const QuaternionFrame& q2) const
{
  QuaternionFrame qi[4] = { QuaternionFrame(q.w(), q.z(), -q.y(), -q.x()),
                            QuaternionFrame(-q.z(), q.w(), q.x(), -q.y()),
                            QuaternionFrame(q.y(), -q.x(), q.w(), -q.z()),
                            QuaternionFrame(q.x(), q.y(), q.z(), q.w()) };
  double dp[4] = { q2.dot(qi[0]), q2.dot(qi[1]), q2.dot(qi[2]), q2.dot(qi[3]) };
  double a[4] = {
    std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]), std::abs(dp[3])
  };

  // find M, N, max indices of a[i]
  int M = 0, N = 1, m = 2, n = 3;
  if (a[M] < a[m])
    std::swap(M, m);
  if (a[N] < a[n])
    std::swap(N, n);
  if (a[M] < a[N]) {
    std::swap(M, N);
    m = n;
  }
  if (a[N] < a[m])
    N = m;

  const double s = std::sqrt(.5f);
  const double h = .5f;

  double vA = a[M];
  double vB = (a[M] + a[N]) * s;
  double vC = (a[0] + a[1] + a[2] + a[3]) * h;

  if (vA > vB && vA > vC) {
    return qi[M].multi(std::copysign(1.f, dp[M]));
  } else if (vB > vC) {
    return qi[M].multi(std::copysign(s, dp[M])) +
           qi[N].multi(std::copysign(s, dp[N]));
  } else {
    return qi[0].multi(std::copysign(h, dp[0])) +
           qi[1].multi(std::copysign(h, dp[1])) +
           qi[2].multi(std::copysign(h, dp[2])) +
           qi[3].multi(std::copysign(h, dp[3]));
  }
}
double
QuaternionFrame::evaluate_closeness(const QuaternionFrame& q1,
                                    const QuaternionFrame& q2)
{
  return Frame3d::metric(q1.to_rotation_matrix(), q2.to_rotation_matrix());
}
}