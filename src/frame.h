#pragma once
#include "eigenDefs.h"

namespace LatticeCore {
class Permutation_3
{
private:
public:
  const std::array<int, 3> data_;

  Permutation_3(int a, int b, int c)
    : data_({ a, b, c }){};

  const int& operator[](int index) const { return data_[index]; }

  /**
   *@brief lv->rv
   */
  Permutation_3 operator*(const Permutation_3& rv) const;

  bool operator==(const Permutation_3& rv) const;

  static const Permutation_3 permutations[6];
  static int transform(int lv, int rv);
};

/**
 *@brief 6-order permutation group
 */
class S6
{
private:
  void constructor(int a, int b, int c);

public:
  S6(int index);
  S6(int a, int b, int c);

  static const S6 TRANS_MATRICES[6];

  Eigen::Matrix<int, 3, 3> trans_matrix_;

  Eigen::Matrix<int, 3, 1> operator*(Eigen::Matrix<int, 3, 1>& rv);

  /**
   *@brief get index z produced by u->v, that M(z) = M(v)*M(u)
   */
  static int index_transform(int u, int v);
};

/**
 * rotation matrix represented frame
 */
class Frame3d
{
private:
  /**
   * u,v must be unit vector, otherwise the result will be wrong
   */
  template<typename T>
  double sin(T u, T v) const
  {
    return (u.cross(v)).norm();
  }

  /**
   *@brief calculate the diffirence after *this do permutation
   */
  double diff(const Frame3d& b, const Permutation_3& permute) const;

public:
  Matrix_3 m_;

  Frame3d()
    : m_(Matrix_3::Identity())
  {}

  Frame3d(const Matrix_3& m)
    : m_(m)
  {}

  Frame3d operator*(const Frame3d& b) const
  {
    return Frame3d(Matrix_3(m_ * b.m_));
  }

  Frame3d operator+(const Frame3d& b) const
  {
    return Frame3d(Matrix_3(m_ + b.m_));
  }

  Matrix_3 to_matrix() const { return m_; }

  /**
   *@brief find the best S6 match between two frames
   */
  int get_S6_matching_index(const Frame3d& b) const;

  double diff(const Frame3d& B) const { return metric(m_, B.m_); }

  // int get_matching_index(const Frame3d &b) { return get_S6_matching_index(b);
  // }

  void normalize();

  /**
    from paper "All-Hex Meshing using Singularity-Restricted Field"
    Yufei Li
  */
  static double harmonic(double x, double y, double z);
  static double harmonic(const V3d& x);
  static double H(const Matrix_3& Rst);
  static Matrix_3 der_H(const Matrix_3& Rst);
  // Partial derivative of function harmonic on x
  static double harmonic_x(double x, double y, double z);
  static double harmonic_x(const V3d& u);
  // Partial derivative of function harmonic on y
  static double harmonic_y(double x, double y, double z);
  static double harmonic_y(const V3d& u);
  // Partial derivative of function harmonic on z
  static double harmonic_z(double x, double y, double z);
  static double harmonic_z(const V3d& u);
  /**
   *@brief distance between two frames A and B.
   */
  static double metric(const Matrix_3& A, const Matrix_3& B);

  static Frame3d average(const Frame3d& a, const Frame3d& b);

  /*********** friend ***************/

  friend std::ostream& operator<<(std::ostream& out, const Frame3d& frame);
};

/**
 * frame class represented by euler angle
 * used in stress field smoothness
 */
class EulerAngleFrame
{
private:
public:
  double u_, v_, w_;

  EulerAngleFrame();
  EulerAngleFrame(double u, double v, double w);
  EulerAngleFrame(const Matrix_3& r);

  Matrix_3 get_matrix();

  /**
   * @brief get rotation matrix from eular angle
   */
  static Matrix_3 R(double u, double v, double w);
  static Matrix_3 dR_1(double u, double v, double w);
  static Matrix_3 dR_2(double u, double v, double w);
  static Matrix_3 dR_3(double u, double v, double w);
};

/**
 * quaternion represented frame
 */
class QuaternionFrame
{
private:
public:
  /**
   * constructor: q(w,x,y,z)
   * coeffients: [x,y,z,w]
   */
  Eigen::Quaterniond q;
  QuaternionFrame()
    : q(Eigen::Quaterniond(1, 0, 0, 0))
  {}
  QuaternionFrame(const Eigen::Quaterniond& p)
    : q(p)
  {}
  QuaternionFrame(const Matrix_3& m)
    : q(m)
  {}
  QuaternionFrame(double x, double y, double z, double w)
    : q(w, x, y, z)
  {}

  double w() const { return q.w(); }
  double x() const { return q.x(); }
  double y() const { return q.y(); }
  double z() const { return q.z(); }

  double& operator[](int i) { return q.coeffs()(i); }

  QuaternionFrame multi(double k) const
  {
    return QuaternionFrame(x() * k, y() * k, z() * k, w() * k);
  }

  QuaternionFrame operator+(const QuaternionFrame& q)
  {
    return QuaternionFrame(x() + q.x(), y() + q.y(), z() + q.z(), w() + q.w());
  }

  QuaternionFrame operator*(const QuaternionFrame& qf2) const
  {
    return QuaternionFrame(q * qf2.q);
  }

  /**
   * @brief get quaternion's coefficients
   *
   * @return Eigen::Vector4d ordered in [x,y,z,w]
   */
  Eigen::Vector4d coeffs() const { return q.coeffs(); }

  QuaternionFrame complement() const;

  Matrix_3 to_rotation_matrix() const;

  double dot(const QuaternionFrame& qf2) const { return q.dot(qf2.q); }

  void normalize() { q.normalize(); }

  void print() const;

  static QuaternionFrame random_quaternion()
  {
   // return QuaternionFrame(Eigen::Quaterniond::UnitRandom());
  }

  static double inner_product(const QuaternionFrame& qf1,
                              const QuaternionFrame& qf2)
  {
    return qf1.q.dot(qf2.q);
  }

  /**
   *@brief find best match of 2 quaternion
   * best match is a quaternion r that makes q1 and q2*r the closest
   */
  static QuaternionFrame request_best_match(const QuaternionFrame& q1,
                                            const QuaternionFrame& q2);

  /**
   *@brief Find the average of two quaternion-represented frames
   *
   * b will be rotated to a. that means quaternion average (a,v*r) will be
   *returned, where r is the best match of q1,q2
   */
  static QuaternionFrame average(const QuaternionFrame& a,
                                 const QuaternionFrame& b);

  /**
   * @brief weighted average of quaternion frames
   */
  static QuaternionFrame weighted_average(
    const std::vector<QuaternionFrame>& list,
    const std::vector<double>& weight,
    bool rotate_invariant = true);

  /**
   * @brief find q1's symmetric rotation that best match q2
   *
   * @param q2
   * @return QuaternionFrame
   */
  QuaternionFrame rotate_to(const QuaternionFrame& q2) const;

  static double evaluate_closeness(const QuaternionFrame& q1,
                                   const QuaternionFrame& q2);
};
}