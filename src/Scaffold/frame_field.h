#pragma once
/**
 * @file frame_field.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief class Permutation_3, S6, Frame3d, EulerAngleFrame, QuaternionFrame
 * @version 0.1
 * @date 2019-10-24
 *
 * @copyright Copyright (c) 2019
 *
 */

#include "kd_tree.hpp"
//#include "eigenDefs.h"
#include "frame.h"

namespace LatticeCore {

class FrameField
{
private:
  std::vector<Frame3d> frames_;

  mutable std::vector<QuaternionFrame>* qframes_ = nullptr;

  std::vector<KDPoint3d> points_;
  std::unique_ptr<KD3d> samples = nullptr;

  void cal_quaternion_frame(std::vector<QuaternionFrame>* f) const;

public:
  FrameField(int n)
    : frames_(n)
  {}

  ~FrameField()
  {
    if (qframes_)
      delete qframes_;
  }

  const Frame3d& operator[](int index) const { return frames_[index]; }
  void resize(int n) { frames_.resize(n); }
  int size() const { return (int)frames_.size(); }
  std::vector<Frame3d>& get_frames() { return frames_; }
  const std::vector<Frame3d>& get_frames() const { return frames_; }

  void set_frames(const std::vector<Matrix_3>& frames);

  std::vector<EulerAngleFrame> request_EulerAngleFrame() const;

  std::vector<QuaternionFrame> request_QuaternionFrame() const;

  const std::vector<QuaternionFrame>& get_qframes() const;

  void correct_frame_direction(int start_cell_id,
                               const SMatrix<int>& adjcent_matrix);

  void build_kd_tree(const std::vector<V3d>& sample_locations);
  void build_kd_tree(const ColMatrix<double>& sample_locations);

  QuaternionFrame interpolate_qframe(const V3d& point) const;
  Frame3d interpolate(const V3d& point) const;
};

} // namespace LatticeCore
