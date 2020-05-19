#include "scaffold.h"
//#include "graph.hpp"
#include "kd_tree.hpp"
#include "log.h"
#include "merge_find_set.hpp"

const static double orthogonal_cosine = std::cos(PI / 12);

namespace LatticeCore {
/***************** local functions *******************/
inline Eigen::Matrix<int, 6, 1>
get_orthogonal_direction_bit(const V3d& rod,
                             const QuaternionFrame& q,
                             double eps)
{
  Matrix_3 Q = q.to_rotation_matrix();
  V3d p = Q.transpose() * rod.normalized();
  Eigen::Matrix<int, 6, 1> res;
  res.setZero();
  for (int j = 0; j < 3; ++j) {
    if ((1.0 - std::abs(p(j))) < eps) {
      if (p(j) > 0) {
        res(j) = 1;
      } else {
        res(j + 3) = 1;
      }
    }
  }
  return res;
}

inline int
get_direction_code(const Eigen::Matrix<int, 6, 1>& bits)
{
  if ((bits.cwiseAbs()).sum() != 1)
    return -1;
  for (int j = 0; j < 6; ++j) {
    if (bits(j) == 1)
      return j;
  }
  return -1;
}

inline bool
most_orthogonal_vertex_valid(const std::vector<int>& candidate_vertex_id,
                             const std::vector<V3d>& candidate_vertex_posision,
                             const V3d& origin,
                             const V3d& ref,
                             int* result)
{
  double max_cos = -1.0;
  V3d r = ref.normalized();
  for (int j : candidate_vertex_id) {
    V3d p = (candidate_vertex_posision[j] - origin).normalized();
    double pcos = p.dot(r);
    if (pcos > max_cos) {
      max_cos = pcos;
      *result = j;
    }
  }
  // too far from orthogonal direction
  if (max_cos < orthogonal_cosine)
    return false;
  else
    return true;
}

inline void
merge_clustered_nodes(const MergeFindSet& findset,
                      const std::vector<V3d>& vertices,
                      const std::vector<IndexPair>& edges,
                      std::vector<V3d>* r_vertices,
                      std::vector<IndexPair>* r_edges)
{
  std::map<int, std::vector<int>> clusters;
  for (int i = 0; i < findset.size(); ++i) {
    clusters[findset.find(i)].push_back(i);
  }
  /* if cluster's size is equal to vertices', then merge is trivial */
  if (clusters.size() == vertices.size()) {
    *r_vertices = vertices;
    *r_edges = edges;
    FLOGOUT(<< "trivial merge")
    return;
  }
  // 表示代表元对应的合并后的点的编号
  std::map<int, int> present_to_clustered;
  // build merged vertices
  int n_merged = clusters.size();
  int k = 0;
  r_vertices->resize(n_merged);
  for (auto& itm : clusters) {
    V3d noval_p(0.0, 0.0, 0.0);
    int min_id = itm.second[0];
    for (auto& id : itm.second) {
      noval_p = noval_p + vertices[id];
      if (id < min_id)
        min_id = id;
    }
    if (itm.second.size() > 1)
      FLOGOUT(<< "min id: " << min_id)
    noval_p /= itm.second.size();
    r_vertices->at(k) = noval_p;
    present_to_clustered[itm.first] = k;
    k++;
  }
  // build new edges
  std::set<IndexPair> merged_edges_set;
  for (auto& e : edges) {
    int u = e[0];
    int v = e[1];
    if (findset.find(u) != findset.find(v)) {
      int noval_u = present_to_clustered[findset.find(u)];
      int noval_v = present_to_clustered[findset.find(v)];
      if (noval_u < noval_v) {
        merged_edges_set.insert(IndexPair{ noval_u, noval_v });
      } else {
        merged_edges_set.insert(IndexPair{ noval_v, noval_u });
      }
    }
  }
  r_edges->assign(merged_edges_set.begin(), merged_edges_set.end());
}

/************ Scaffold **************************/

void
Scaffold::remove_redundancy(const std::vector<V3d>& points,
                            const std::vector<IndexPair>& edges,
                            std::vector<V3d>* out_points,
                            std::vector<IndexPair>* out_edges)
{
  FLOGOUT(<< "remove redundancy -> initial vertices,edges " << points.size()
          << ',' << edges.size())
  /* stage 1 */
  std::vector<IndexPair> remaining_edges;
  int n_points = points.size();
  MergeFindSet merge_set(n_points);
  for (auto& e : edges) {
    int u = e[0];
    int v = e[1];
    double dist = (points[u] - points[v]).norm();
    if (dist < this->vertex_merge_threshold) {
      merge_set.merge(u, v);
    } else {
      remaining_edges.push_back({ u, v });
    }
  }
  std::vector<V3d> stage1_points;
  std::vector<IndexPair> stage1_edges;
  merge_clustered_nodes(
    merge_set, points, remaining_edges, &stage1_points, &stage1_edges);
  FLOGOUT(<< "remove redundancy -> stage 1 vertices,edges "
          << stage1_points.size() << ',' << stage1_edges.size())
  /* stage 2 */
  /* merge vertices that are close but have no common edge */
  int n_points2 = stage1_points.size();
  MergeFindSet merge_set2(n_points2);
  KD3d kdtree(stage1_points);
  for (int i = 0; i < n_points2; ++i) {
    auto d =
      kdtree.radiusSearch(stage1_points[i], this->vertex_merge_threshold);
    for (int j : d) {
      merge_set2.merge(i, j);
    }
  }
  merge_clustered_nodes(
    merge_set2, stage1_points, stage1_edges, out_points, out_edges);
  FLOGOUT(<< "remove redundancy -> stage 2 vertices,edges "
          << out_points->size() << ',' << out_edges->size())
}

void
Scaffold::raw_graph_merge()
{
  remove_redundancy(
    *raw_vertices, *raw_edges, &scaffold_vertices, &scaffold_edges);
}

void
Scaffold::repair_vertices()
{
  cal_merged_qframes();
  HierarchyLayer dense_layer;
  dense_layer.vertices = scaffold_vertices;
  dense_layer.qframes = scaffold_qframes;
  dense_layer.edge_list = scaffold_edges;
  dense_layer.increase_vertex_density(ff, cell_scale);
  auto optimized_positions = dense_layer.vertices;
#ifdef LATTICE_USE_TBB
  int iteration_time = 300;
  dense_layer.parallel_optimize_positions(
    &optimized_positions, iteration_time, cell_scale);
#else
  int iteration_time = 70;
  dense_layer.optimize_positions(
    &optimized_positions, iteration_time, cell_scale);
#endif
  remove_redundancy(optimized_positions,
                    dense_layer.edge_list,
                    &scaffold_vertices,
                    &scaffold_edges);
  cal_merged_qframes();
}

std::set<int>
Scaffold::insert_new_vertex(const KD3d& old_vertices,
                            const std::vector<V3d>& candidate_vertices,
                            const V3d& new_vertex)
{
  double insert_threshold = this->cell_scale * 0.3;
  std::set<int> new_edges;
  int old_vertices_num = (int)scaffold_vertices.size();
  auto adj_vertex_ids =
    old_vertices.radiusSearch(new_vertex, edge_connect_threshold);
  // existing vertices are close to new_vertex, then cancel insert
  for (int j : adj_vertex_ids) {
    double dist = (scaffold_vertices[j] - new_vertex).norm();
    if (dist < insert_threshold) {
      return std::set<int>();
    }
  }
  std::vector<int> candidate_in_neighborhood;
  for (int i = 0; i < candidate_vertices.size(); ++i) {
    double dist = (candidate_vertices[i] - new_vertex).norm();
    if (dist < vertex_merge_threshold) {
      return std::set<int>();
    } else if (dist < edge_connect_threshold) {
      candidate_in_neighborhood.push_back(i);
    }
  }

  // connect orthogonal vertices
  int closest_vertex_id;
  auto Q = ff->interpolate(new_vertex).to_matrix();
  for (int i = 0; i < 3; ++i) {
    V3d ref = Q.col(i);
    if (most_orthogonal_vertex_valid(adj_vertex_ids,
                                     this->scaffold_vertices,
                                     new_vertex,
                                     ref,
                                     &closest_vertex_id)) {
      new_edges.insert(closest_vertex_id);
    }
    if (most_orthogonal_vertex_valid(candidate_in_neighborhood,
                                     candidate_vertices,
                                     new_vertex,
                                     ref,
                                     &closest_vertex_id)) {
      new_edges.insert(closest_vertex_id + old_vertices_num);
    }
    ref = -Q.col(i);
    if (most_orthogonal_vertex_valid(adj_vertex_ids,
                                     this->scaffold_vertices,
                                     new_vertex,
                                     ref,
                                     &closest_vertex_id)) {
      new_edges.insert(closest_vertex_id);
    }
    if (most_orthogonal_vertex_valid(candidate_in_neighborhood,
                                     candidate_vertices,
                                     new_vertex,
                                     ref,
                                     &closest_vertex_id)) {
      new_edges.insert(closest_vertex_id + old_vertices_num);
    }
  }
  return new_edges;
}

void
Scaffold::repair_connection()
{
  int n_vertices = scaffold_vertices.size();
  kd_ptr = std::make_unique<KD3d>(scaffold_vertices);
  KD3d& kdtree = *kd_ptr;
  auto adjset = get_adj_list();

  double radius = cell_scale * config_table.connect_rate;
  double eps = std::cos(PI / 6);
  std::set<IndexPair> complement_edges;
  std::vector<V3d> added_vertex_position;
  std::vector<std::set<int>> added_vertex_adj_list;

  auto complete_vertex_neighbor =
    [&](std::vector<int>& d, int i, const V3d& nor) {
      // find an orthogonal vertex that is not connected
      V3d origin = scaffold_vertices[i];
      // FLOGOUT(<<"nor: "<<nor.transpose())
      int candidate_id;
      // if candidate point is an existing node
      if (most_orthogonal_vertex_valid(
            d, this->scaffold_vertices, origin, nor, &candidate_id) &&
          (adjset[i].find(candidate_id) == adjset[i].end())) {
        if (candidate_id > i) {
          complement_edges.insert({ i, candidate_id });
        } else {
          complement_edges.insert({ candidate_id, i });
        }
      }
    };

  patch_vertices.clear();
  for (int i = 0; i < n_vertices; ++i) {
    auto d = kdtree.radiusSearch(scaffold_vertices[i], radius);
    // in neighborhood but have not an edge
    Matrix_3 Q = scaffold_qframes[i].to_rotation_matrix();
    for (int k = 0; k < 3; ++k) {
      complete_vertex_neighbor(d, i, Q.col(k) * cell_scale);
      complete_vertex_neighbor(d, i, -Q.col(k) * cell_scale);
    }
  }
  this->patch_vertices = added_vertex_position;
  scaffold_edges.insert(
    scaffold_edges.end(), complement_edges.begin(), complement_edges.end());
  FLOGOUT(<< "edges after repair: " << scaffold_edges.size())
}

void
Scaffold::remove_bevel_edges()
{
  std::vector<IndexPair> remaining_edges;
  for (auto& e : scaffold_edges) {
    int u = e[0];
    int v = e[1];
    V3d rod = scaffold_vertices[u] - scaffold_vertices[v];
    double dist = rod.norm();
    if (dist < edge_connect_threshold) {
      const auto& q1 = scaffold_qframes[u];
      QuaternionFrame q2 = scaffold_qframes[v].rotate_to(q1);
      Matrix_3 Q = QuaternionFrame::average(q1, q2).to_rotation_matrix();
      V3d t = Q.transpose() * rod;
      int none_zero = 3;
      for (int j = 0; j < 3; ++j) {
        if (std::abs(t(j)) < vertex_merge_threshold)
          none_zero--;
      }
      if (none_zero == 1) {
        remaining_edges.push_back({ u, v });
      }
    }
  }
  scaffold_edges = std::move(remaining_edges);
}

void
Scaffold::cal_vertex_irregularity()
{
  int n_vertices = scaffold_vertices.size();
  get_adj_list();
  double eps = std::sin(PI / 8);
  vertex_regularity.resize(n_vertices);
  orthogonal_neighbor.resize(n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    vertex_regularity[i] = 0.0;
    Eigen::Matrix<int, 6, 1> neighbor_bit;
    neighbor_bit.setZero();
    for (int j : adj_list[i]) {
      V3d rod = scaffold_vertices[i] - scaffold_vertices[j];
      auto bit = get_orthogonal_direction_bit(rod, scaffold_qframes[i], eps);
      neighbor_bit += bit;
      // if (bit.sum() == 0) {
      // vertex_regularity[i] += 1.0;
      // }
    }
    for (int j = 0; j < 6; ++j) {
      if (neighbor_bit(j) > 1 || neighbor_bit(j) == 0) {
        vertex_regularity[i] += 1;
      }
    }
    orthogonal_neighbor[i] = neighbor_bit;
  }
}

const std::vector<std::set<int>>&
Scaffold::get_adj_list() const
{
  int n_vertices = scaffold_vertices.size();
  if (adj_list.size() < n_vertices) {
    adj_list = edge_pair_to_adjacent_vector(scaffold_edges, n_vertices);
  }
  return adj_list;
}

void
Scaffold::cal_merged_qframes()
{
  int n = scaffold_vertices.size();
  scaffold_qframes.clear();
  scaffold_qframes.reserve(n);
  for (int i = 0; i < n; ++i) {
    scaffold_qframes.push_back(ff->interpolate_qframe(scaffold_vertices[i]));
  }
}

void
Scaffold::cluster_regular_node()
{
  int n = scaffold_vertices.size();
  std::vector<bool> visited(n, false);
  std::queue<int> list;
  for (int ka = 0; ka < n; ++ka) {
    if (visited[ka])
      continue;
    list.push(ka);
    visited[ka] = true;
    get_adj_list();
    std::vector<int> bundle;
    while (!list.empty()) {
      int u = list.front();
      list.pop();
      if (vertex_regularity[u] > 0.0) {
        continue;
      }
      bundle.push_back(u);
      for (int v : adj_list[u]) {
        if (!visited[v]) {
          list.push(v);
          visited[v] = true;
        }
      }
    }
    this->clusters.push_back(bundle);
  }
}

Scaffold::Scaffold(Hierarchy* pyramid, FrameField* field)
  : pyr(pyramid)
  , ff(field)
{
  auto& layer = pyramid->layers[0];
  // field->build_kd_tree(layer.vertices);
  raw_vertices = &(pyramid->optimized_positions);
  raw_edges = &(layer.edge_list);
  raw_frames = &(layer.qframes);
  cell_scale = pyramid->get_cell_scale();
  FLOGOUT(<< "scaffold cell scale: " << cell_scale)
  int num_edges = layer.edge_list.size();
  vertex_merge_threshold = cell_scale * config_table.merge_rate;
  edge_connect_threshold = cell_scale * config_table.connect_rate;
}

void
Scaffold::postprocess()
{
  raw_graph_merge();
  repair_vertices();
  // if (scaffold_qframes.size() != scaffold_vertices.size()) {
  //   cal_merged_qframes();
  // }
  cal_vertex_irregularity();
  repair_connection();
  remove_bevel_edges();
  // cluster_regular_node();
}
}