#pragma once
#include "frame_field.h"
#include "eigenDefs.h"

namespace LatticeCore {
class Scaffold
{
private:
  const FrameField* ff = nullptr;
  //const Hierarchy* pyr = nullptr;
  const std::vector<V3d>* raw_vertices;
  const std::vector<IndexPair>* raw_edges;
  const std::vector<QuaternionFrame>* raw_frames;
  mutable std::vector<std::set<int>> adj_list;
  std::vector<Eigen::Matrix<int, 6, 1>> orthogonal_neighbor;
  double vertex_merge_threshold;
  double edge_connect_threshold;
  std::unique_ptr<KD3d> kd_ptr = nullptr;

  struct ScaffoldConfig
  {
    double merge_rate = 0.2;
    double connect_rate = 1.21;
  } config_table;

  const std::vector<std::set<int>>& get_adj_list() const;

  void cal_merged_qframes();

  void cluster_regular_node();

  void remove_redundancy(const std::vector<V3d>& points,
                         const std::vector<IndexPair>& edges,
                         std::vector<V3d>* out_points,
                         std::vector<IndexPair>* out_edges);

  void raw_graph_merge();

  void repair_vertices();

  std::set<int> insert_new_vertex(const KD3d& kdtree,
                                  const std::vector<V3d>& candidate_vertices,
                                  const V3d& new_vertex);

  /**
   * @brief add vertices in proper positions to get a better mesh structure
   */
  void repair_connection();

  /**
   */
  void remove_bevel_edges();

  void cal_vertex_irregularity();

public:
  std::vector<V3d> scaffold_vertices;
  std::vector<QuaternionFrame> scaffold_qframes;
  std::vector<IndexPair> scaffold_edges;
  std::vector<double> vertex_regularity;
  std::vector<V3d> patch_vertices; // points added to generate hexahedra
  std::vector<std::vector<int>> clusters;
  double cell_scale;

 // Scaffold(Hierarchy* pyramid, FrameField* field);

  void set_merge_threshold(double threshold)
  {
    config_table.merge_rate = threshold;
    vertex_merge_threshold = cell_scale * threshold;
  }

  void postprocess();
};
}