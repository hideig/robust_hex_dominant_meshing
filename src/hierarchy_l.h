#pragma once
/**
 * @file hierarchy.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief
 * @version 0.1
 * @date 2020-02-23
 *
 * @copyright Copyright (c) 2020
 *
 */

#include "frame_field.h"
#include "dual_graph.hpp"
#include "eigenDefs.h"

namespace LatticeCore {

class HierarchyLayer
{
private:
public:
  /**************** data ***************/

  std::vector<V3d> vertices;
  std::vector<QuaternionFrame> qframes;

  /**
   * @brief adjacent matrix
   *
   * if link_matrix(i,j) != 0, then i,j have common edge
   * link_matrix(i,j) is the difference of i,j's frames + 1,
   * +1 is to avoid 0 term disappear when i,j's frames are identical
   */
  // SMatrix<double> link_matrix;
  std::vector<IndexPair> edge_list;

  /*********** functions *****************/

  void set_vertices(const std::vector<V3d>& vertices);

  inline double edge_diff(const IndexPair& e)
  {
    return QuaternionFrame::evaluate_closeness(qframes[e[0]], (qframes[e[1]]));
  }

  void sort_edges_by_qframe_diff();

  std::vector<V3d> increase_vertex_density(const FrameField* ff, double radius);

  void optimize_positions(std::vector<V3d>* node_position_ptr,
                          int iteration_times,
                          double cell_scale);

#ifdef LATTICE_USE_TBB
  void parallel_optimize_positions(std::vector<V3d>* node_position_ptr,
                                   int iteration_times,
                                   double cell_scale);
#endif

  std::vector<IndexPair> get_edges() const;
};

class Hierarchy
{
private:
  double cell_scale;
  HierarchyLayer dense_layer;

  /********** functions **************/

  std::vector<std::pair<int, int>> collapse_vertices(
    HierarchyLayer& lower_layer,
    std::vector<V3d>* vertices,
    std::vector<QuaternionFrame>* qframe);

  std::vector<IndexPair> build_higher_edges(
    const std::vector<std::pair<int, int>>& remaining_edges,
    const std::vector<QuaternionFrame>& higher_frames);

  bool build_layer();

  /**
   * @brief cast layer_id's positions to lower layer
   *
   * @param layer_id id of upper layer
   */
  void positions_down_cast(std::vector<V3d>* positions,

                           const std::vector<int>& cast_table);

  void extract_tau();

public:
  std::vector<HierarchyLayer> layers;
  std::vector<V3d> optimized_positions;
  std::vector<V3d> noval_vertices;

  // if vertices is less than min_vertices,
  // then hierarchy layer construction will stop
  int min_vertices = 10;

  // std::vector<SMatrix<double>> upcast_matrix;
  std::vector<std::vector<int>> upcast_tables;
  std::vector<Eigen::Vector3d> taus;

  int build_hierarchy(const DualGraph& dg,
                      const FrameField& ff,
                      double scale,
                      int depth);

  double get_cell_scale() { return cell_scale; }
  void set_cell_scale(double cell_size)
  {
    if (cell_size > 0.0)
      this->cell_scale = cell_size;
  }

  void optimize();

#ifdef LATTICE_USE_TBB
  void parallel_optimize();
#endif

  /**
   * add more vertices to get a more regular mesh
   */
  void dense_optimize(const FrameField* ff);
};
} // namespace LatticeCore
