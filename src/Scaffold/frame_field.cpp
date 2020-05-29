#include "frame_field.h"
#include <cmath>
#include <iostream>


namespace LatticeCore {

void
FrameField::cal_quaternion_frame(std::vector<QuaternionFrame>* f) const
{
  int n = frames_.size();
  f->clear();
  f->reserve(n);
  for (auto& m : frames_) {
    f->push_back(QuaternionFrame(m.to_matrix()));
  }
}

void
FrameField::set_frames(const std::vector<Matrix_3>& new_frames)
{
  int n = new_frames.size();
  frames_.resize(n);
  for (int i = 0; i < n; ++i) {
    frames_[i] = new_frames[i];
  }
  if (!qframes_) {
    qframes_ = new std::vector<QuaternionFrame>;
  }
  cal_quaternion_frame(qframes_);
}

std::vector<EulerAngleFrame>
FrameField::request_EulerAngleFrame() const
{
  std::vector<EulerAngleFrame> frames;
  frames.reserve(size());

  for (auto& c : frames_) {
    frames.push_back(EulerAngleFrame(c.m_));
  }
  return frames;
}

std::vector<QuaternionFrame>
FrameField::request_QuaternionFrame() const
{
  std::vector<QuaternionFrame> f;
  cal_quaternion_frame(&f);
  return f;
}

const std::vector<QuaternionFrame>&
FrameField::get_qframes() const
{
  if (!qframes_) {
    qframes_ = new std::vector<QuaternionFrame>;
    cal_quaternion_frame(qframes_);
  }
  return *qframes_;
}

void
FrameField::correct_frame_direction(int start_cell_id,
                                    const SMatrix<int>& adjcent_matrix)
{
  // timer
  //DebugTools::Log::get_instance().tick();

  // bfs correct frames' major directions
  int n = size();
  std::vector<bool> visited(n, false);
  std::queue<int> bfs_list;
  // bfs_list.push(start_cell_id);
  // visited[start_cell_id] = true;
  for (int cell_id = 0; cell_id < n; ++cell_id) {
    if (!visited[cell_id]) {
      bfs_list.push(cell_id);
      visited[cell_id] = true;
      while (!bfs_list.empty()) {
        int u = bfs_list.front();
        bfs_list.pop();
        for (SMatrix<int>::InnerIterator it(adjcent_matrix, u); it; ++it) {
          int v = it.col();
          if (!visited[v]) {
            V3d u_major_vector = frames_[u].m_.col(0);
            V3d v_major_vector = frames_[v].m_.col(0);
            double cosin = u_major_vector.dot(v_major_vector);
            if (std::abs(cosin) > 1 - 1e-1) {
              if (cosin < 0)
                frames_[v].m_ *= -1;
              bfs_list.push(v);
              visited[v] = true;
            }
          } else {
            V3d u_major_vector = frames_[u].m_.col(0);
            V3d v_major_vector = frames_[v].m_.col(0);
            double cosin = u_major_vector.dot(v_major_vector);
            if (std::abs(cosin + 1) < 1e-2) {
             // FLOGOUT(<< "direction cofliction! " << u << ' ' << v)
            }
          }
        }
      }
    }
  }
  int cnt = 0;
  for (bool j : visited) {
    if (j)
      cnt++;
  }
 // FLOGOUT(<< "total corrected frames: " << cnt << " totoal frames: " << n
    //      << std::endl)
  // timer
 // DebugTools::Log::get_instance().tick();
 /* std::cout << "correction time cosumed: "
            << DebugTools::Log::get_instance().get_millisec_span()
            << " millisec" << std::endl;*/
}

void
FrameField::build_kd_tree(const std::vector<V3d>& sample_locations)
{
  int n = sample_locations.size();
  if (n != frames_.size()) {
    std::cerr << "sample points' number is not equal to frames' number!"
              << std::endl;
  }
  samples = std::make_unique<KD3d>(sample_locations);
}

void
FrameField::build_kd_tree(const ColMatrix<double>& sample_locations)
{
  int n = sample_locations.cols();
  if (n != frames_.size()) {
    std::cerr << "sample points' number is not equal to frames' number!"
              << std::endl;
  }
  samples = std::make_unique<KD3d>(sample_locations);
}

QuaternionFrame
FrameField::interpolate_qframe(const V3d& point) const
{
  if (!samples) {
    std::cerr << "kdtree has not been built!" << std::endl;
    return QuaternionFrame();
  }
  int research_num = 6;
  auto result = samples->knnSearch(point, research_num);
  auto& qframes = get_qframes();

  // interpolate frame from neighboring research_num samples
  std::vector<QuaternionFrame> qf;
  std::vector<double> weight;
  for (int id : result) {
    double dist = (point - samples->get_point(id)).norm();
    qf.push_back(qframes[id]);
    weight.push_back(std::exp(-dist * dist));
  }
  QuaternionFrame r = QuaternionFrame::weighted_average(qf, weight);
  return r;
}

Frame3d
FrameField::interpolate(const V3d& point) const
{
  return interpolate_qframe(point).to_rotation_matrix();
}
} // namespace LatticeCore