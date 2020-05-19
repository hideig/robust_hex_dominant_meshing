#include "hierarchy.h"
#include "LatticeCore/Tools/merge_find_set.hpp"
#include <algorithm>
#include <random>
#ifdef LATTICE_USE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#endif

// for debug
#include "LatticeCore/Tools/log.h"

namespace LatticeCore {

inline V3d
cal_tau(const V3d& p1,
        const V3d& p_ref,
        const QuaternionFrame& q1,
        const QuaternionFrame& q_ref)
{
  Matrix_3 Q = QuaternionFrame::average(q_ref, q1).to_rotation_matrix();
  return Q.transpose() * (p1 - p_ref);
}

inline bool
is_orthogonal_vector(const V3d& tau)
{
  Eigen::Vector3d t = Eigen::round(tau.array());
  int ones = 0;
  for (int j = 0; j < 3; ++j) {
    int i = static_cast<int>(t(j) < 0 ? -t(j) : t(j));
    if (i == 1)
      ones++;
    else if (i > 0)
      return false;
  }
  return ones == 1;
}

/**
 * 求离p最近的ref_p的整数点
 */
inline V3d
integer_move(V3d p,
             V3d ref_p,
             QuaternionFrame& q,
             QuaternionFrame& ref_q,
             double scale,
             double invScale)
{
  auto q_j = q.rotate_to(ref_q);
  Matrix_3 Q = QuaternionFrame::average(ref_q, q_j).to_rotation_matrix();
  V3d tau_d = Q.transpose() * (p - ref_p) * invScale;
  tau_d = Eigen::round(tau_d.array());
  return ref_p + Q * tau_d * scale;
}

void
HierarchyLayer::set_vertices(const std::vector<V3d>& vertices)
{
  int n = vertices.size();
  this->vertices.resize(n);
  for (int i = 0; i < n; ++i) {
    this->vertices[i] = vertices[i];
  }
}

void
HierarchyLayer::sort_edges_by_qframe_diff()
{
  std::sort(edge_list.begin(),
            edge_list.end(),
            [this](const IndexPair& e1, const IndexPair& e2) {
              return this->edge_diff(e1) < this->edge_diff(e2);
            });
}

std::vector<V3d>
HierarchyLayer::increase_vertex_density(const FrameField* ff, double radius)
{
  KD3d kdtree(vertices);
  int old_vertices_n = vertices.size();
  std::vector<std::vector<int>> adjlist(old_vertices_n);
  for (auto& e : edge_list) {
    adjlist[e[0]].push_back(e[1]);
    adjlist[e[1]].push_back(e[0]);
  }
  const double merge_criterion = radius * 0.3;
  const double connect_criterion = radius * 1.2;
  std::vector<V3d> noval_vertices;
  std::vector<IndexPair> noval_edges;
  std::vector<std::vector<int>> old_v_nv_neighborhood(old_vertices_n);
  auto insert_new_vertex = [&](int i, const V3d& p) {
    // if is not near an existing vertex, then new vertex is valid
    for (auto& k : adjlist[i]) {
      if ((p - vertices[k]).norm() < merge_criterion)
        return;
    }
    // return valid;
    auto nearby_old_vertices = kdtree.radiusSearch(p, connect_criterion);
    if (nearby_old_vertices.size() < 2)
      return;
    // no near existing new vertices
    for (int k : nearby_old_vertices) {
      for (int j : old_v_nv_neighborhood[k]) {
        if ((p - noval_vertices[j - old_vertices_n]).norm() < merge_criterion)
          return;
      }
    }
    noval_vertices.push_back(p);
    int noval_id = old_vertices_n + noval_vertices.size() - 1;
    for (int k : nearby_old_vertices) {
      noval_edges.push_back({ k, noval_id });
      old_v_nv_neighborhood[k].push_back(noval_id);
    }
  };
  int noval_v_num = 0;
  for (int i = 0; i < old_vertices_n; ++i) {
    Matrix_3 Q = qframes[i].to_rotation_matrix();
    const V3d& origin = vertices[i];
    for (int j = 0; j < 3; ++j) {
      V3d r = origin + Q.col(j) * radius;
      insert_new_vertex(i, r);
      r = origin - Q.col(j) * radius;
      insert_new_vertex(i, r);
    }
  }
  KD3d noval_knn(noval_vertices);
  for (int i = 0; i < noval_vertices.size(); ++i) {
    V3d& p = noval_vertices[i];
    auto d = noval_knn.radiusSearch(p, merge_criterion);
    for (int j : d) {
      if (i < j)
        noval_edges.push_back({ i + old_vertices_n, j + old_vertices_n });
    }
  }
  std::vector<QuaternionFrame> noval_f;
  for (auto& v : noval_vertices) {
    noval_f.push_back(ff->interpolate_qframe(v));
  }

  FLOGOUT(<< "dense layer -> vertices: " << vertices.size() << "  "
          << noval_vertices.size());
  vertices.insert(vertices.end(), noval_vertices.begin(), noval_vertices.end());

  FLOGOUT(<< "dense layer -> qframes: " << qframes.size() << "  "
          << noval_f.size());
  qframes.insert(qframes.end(), noval_f.begin(), noval_f.end());

  FLOGOUT(<< "dense layer -> edges: " << edge_list.size() << "  "
          << noval_edges.size());
  edge_list.insert(edge_list.end(), noval_edges.begin(), noval_edges.end());
  return noval_vertices;
}

void
HierarchyLayer::optimize_positions(std::vector<V3d>* node_position_ptr,
                                   int iteration_times,
                                   double cell_scale)
{
  int n = vertices.size();
  double invScale = 1.0 / cell_scale;
  // ColMatrix<double> base_positions = positions;
  std::vector<V3d>& positions = *node_position_ptr;
  std::random_device rd;
  auto rd_engine = std::default_random_engine{ rd() };
  // build adjacent lists
  std::vector<std::vector<std::pair<int, double>>> neighbors(n);
  for (auto& e : edge_list) {
    int u = e[0];
    int v = e[1];
    double dist = (vertices[u] - vertices[v]).norm();
    if (dist < 1e-2)
      dist = 1e-2;
    double value = cell_scale / dist;
    neighbors[u].push_back({ v, value });
    neighbors[v].push_back({ u, value });
  }

  // optimize every nodes
  std::vector<int> node_ids(n);
  for (int i = 0; i < n; ++i)
    node_ids[i] = i;
  while (iteration_times--) {
    // for (int i = 0; i < n; ++i) {
    std::shuffle(node_ids.begin(), node_ids.end(), rd_engine);
    for (int i : node_ids) {
      const auto& q_i = qframes[i];
      V3d p_i = positions.at(i);
      V3d nov_p = { 0, 0, 0 };
      double total_weight = 0.0;
      // shuffle
      std::shuffle(neighbors[i].begin(), neighbors[i].end(), rd_engine);
      // double k = 0.0;
      for (auto& nei_pair : neighbors[i]) {
        int j = nei_pair.first;
        double weight = nei_pair.second;
        V3d p_j = positions.at(j);
        V3d o_j =
          integer_move(p_i, p_j, qframes[i], qframes[j], cell_scale, invScale);
        nov_p += o_j * weight;
        total_weight += weight;
      }
      p_i = nov_p / total_weight;
      // move nov_p closer to vertex
      Matrix_3 Q = q_i.to_rotation_matrix();
      V3d tau = Q.transpose() * (vertices[i] - p_i) * invScale;
      // V3d tau = Q.transpose() * (base_positions.col(i) - p_i) * invScale;
      tau = Eigen::round(tau.array());
      p_i = p_i + Q * tau * cell_scale;
      // update neighbor weights
      for (auto& nei_pair : neighbors[i]) {
        double dist = (p_i - positions.at(nei_pair.first)).norm();
        if (dist < 1e-2)
          dist = 1e-2;
        nei_pair.second = cell_scale / dist;
      }
      positions.at(i) = p_i;
    }
  }
}

#ifdef LATTICE_USE_TBB
void
HierarchyLayer::parallel_optimize_positions(std::vector<V3d>* node_position_ptr,
                                            int iteration_times,
                                            double cell_scale)
{
  int n_vertices = vertices.size();
  double invScale = 1.0 / cell_scale;
  auto& positions = *node_position_ptr;

  while (iteration_times--) {
    std::vector<V3d> optimized_p(positions.size());
    // build adjacent lists
    std::vector<std::vector<std::pair<int, double>>> neighbors(n_vertices);
    for (auto& e : edge_list) {
      int u = e[0];
      int v = e[1];
      double dist = (positions[u] - positions[v]).norm();
      if (dist < 1e-2)
        dist = 1e-2;
      double value = cell_scale / dist;
      neighbors[u].push_back(std::make_pair(v, value));
      neighbors[v].push_back(std::make_pair(u, value));
    }

    tbb::parallel_for(
      tbb::blocked_range<int>(0, n_vertices, 128),
      [&](const tbb::blocked_range<int>& range) {
        std::random_device rd;
        auto rd_engine = std::default_random_engine{ rd() };
        // optimize every nodes
        for (int i = range.begin(); i != range.end(); ++i) {
          std::shuffle(neighbors[i].begin(), neighbors[i].end(), rd_engine);
          const auto& q_i = qframes[i];
          V3d p_i = positions.at(i);
          V3d nov_p = { 0, 0, 0 };
          double total_weight = 0.0;
          // double k = 0.0;
          for (auto& nei_pair : neighbors[i]) {
            int j = nei_pair.first;
            double weight = nei_pair.second;
            V3d p_j = positions.at(j);
            V3d o_j = integer_move(
              p_i, p_j, qframes[i], qframes[j], cell_scale, invScale);
            nov_p += o_j * weight;
            total_weight += weight;
          }
          p_i = nov_p / total_weight;
          // move nov_p closer to vertex
          Matrix_3 Q = q_i.to_rotation_matrix();
          V3d tau = Q.transpose() * (vertices[i] - p_i) * invScale;
          // V3d tau = Q.transpose() * (base_positions.col(i) - p_i) *
          // invScale;
          tau = Eigen::round(tau.array());
          p_i = p_i + Q * tau * cell_scale;
          optimized_p.at(i) = p_i;
        }
      });

    positions = std::move(optimized_p);
  }
}
#endif

std::vector<IndexPair>
HierarchyLayer::get_edges() const
{
  int n = edge_list.size();
  std::vector<std::array<IndexType, 2>> e(n);
  for (int i = 0; i < n; ++i) {
    e[i][0] = edge_list[i][0];
    e[i][1] = edge_list[i][1];
  }
  return e;
}

std::set<std::tuple<double, int, int>>
matrix_to_set(const SMatrix<double>& adj_M)
{
  std::set<std::tuple<double, int, int>> edge_set;
  edge_set.clear();
  // construct edge_pairs from adjacent matrix
  for (int i = 0; i < adj_M.outerSize(); ++i) {
    for (SMatrix<double>::InnerIterator it(adj_M, i); it; ++it) {
      int row = it.row();
      int col = it.col();
      if (row >= col)
        break;
      // QuaternionFrame &f1 = lower_layer.frames[row];
      // QuaternionFrame &f2 = lower_layer.frames[col];
      edge_set.insert(std::make_tuple(it.value(), row, col));
    }
  }
  return edge_set;
}

/******************* Hierarchy **************************/

std::vector<std::pair<int, int>>
Hierarchy::collapse_vertices(HierarchyLayer& lower_layer,
                             std::vector<V3d>* vertices,
                             std::vector<QuaternionFrame>* qframe)
{
  int lower_n = lower_layer.vertices.size();
  FLOGOUT(<< "lower_n: " << lower_n)
  lower_layer.sort_edges_by_qframe_diff();
  auto& edge_set = lower_layer.edge_list;

  std::vector<bool> edge_collapsed(edge_set.size(), false);
  // higer layer vertice to lower vertices
  // std::vector<Tripled> upcast_vertices_triples;
  std::vector<int> upcast_t(lower_n, -1);
  // higer layer vertice to lower vertices
  std::vector<V3d> higher_vertices;
  auto& higher_frames = *qframe;
  // higher edges are claimed after new edges's triples constructed
  int higher_vertices_n = 0;
  higher_vertices.clear();
  // construct new edge_pairs
  std::vector<std::pair<int, int>> remaining_edges;
  remaining_edges.clear();
  std::vector<bool> collpased(lower_n, false);

  for (auto& e : edge_set) {
    int u = e[1];
    int v = e[0];
    if (!collpased[u] && !collpased[v]) {
      collpased[u] = collpased[v] = true;
      QuaternionFrame frame = QuaternionFrame::average(lower_layer.qframes[u],
                                                       lower_layer.qframes[v]);
      higher_vertices.push_back(
        (lower_layer.vertices[u] + lower_layer.vertices[v]) / 2);
      higher_frames.push_back(frame);
      higher_vertices_n = higher_vertices.size() - 1;
      upcast_t[u] = higher_vertices_n;
      upcast_t[v] = higher_vertices_n;
    } else {
      // not collapsed edges will be the higer layer's edges
      remaining_edges.push_back(std::make_pair(u, v));
    }
  }

  // deal with un-collapsed vertices
  for (int i = 0; i < lower_n; ++i) {
    if (!collpased[i]) {
      higher_vertices.push_back(lower_layer.vertices[i]);
      higher_frames.push_back(lower_layer.qframes[i]);
      higher_vertices_n++;
      upcast_t[i] = higher_vertices_n;
    }
  }
  higher_vertices_n++;
  vertices->resize(higher_vertices_n);
  // store vertices
  for (int i = 0; i < higher_vertices_n; ++i) {
    vertices->at(i) = higher_vertices[i];
  }

  // save up/down cast matrix
  this->upcast_tables.push_back(upcast_t);

  return remaining_edges;
}

std::vector<IndexPair>
Hierarchy::build_higher_edges(
  const std::vector<std::pair<int, int>>& remaining_edges,
  const std::vector<QuaternionFrame>& higher_frames)
{
  auto& upcast_t = upcast_tables[upcast_tables.size() - 1];
  std::set<std::tuple<int, int>> connected_vertices;
  std::vector<IndexPair> higher_edges;
  for (auto& e : remaining_edges) {
    int u = upcast_t[e.first];
    int v = upcast_t[e.second];
    // u always less than v
    if (v < u)
      std::swap(u, v);
    if (connected_vertices.find(std::tie(u, v)) == connected_vertices.end()) {
      higher_edges.push_back({ u, v });
      connected_vertices.insert(std::make_tuple(u, v));
    }
  }
  return higher_edges;
}

bool
Hierarchy::build_layer()
{
  // collapse vertices
  auto& lower_layer = layers[layers.size() - 1];
  auto newlayer = HierarchyLayer();
  auto remaining_edges =
    collapse_vertices(lower_layer, &(newlayer.vertices), &(newlayer.qframes));
  int higher_vertices_n = newlayer.vertices.size();
  FLOGOUT(<< "vertices: " << higher_vertices_n)
  newlayer.edge_list = build_higher_edges(remaining_edges, newlayer.qframes);
  layers.push_back(newlayer);
  if (higher_vertices_n <= this->min_vertices) {
    return false;
  }
  return true;
}

void
Hierarchy::positions_down_cast(std::vector<V3d>* positions,
                               const std::vector<int>& cast_table)
{
  int n = cast_table.size();
  std::vector<V3d> new_p(n);
  for (int i = 0; i < n; ++i) {
    int uid = cast_table[i];
    new_p.at(i) = positions->at(uid);
  }
  *positions = std::move(new_p);
}

void
Hierarchy::extract_tau()
{
  HierarchyLayer& layer = layers[0];
  const auto& p = optimized_positions;
  int num_edges = layer.edge_list.size();
  double invScale = 1.0 / cell_scale;

  // get tau of each edges
  int i = 0;
  taus.resize(num_edges);
  for (auto e : layer.edge_list) {
    int u = e[0];
    int v = e[1];
    QuaternionFrame q1 = layer.qframes[u];
    QuaternionFrame q2 = layer.qframes[v].rotate_to(layer.qframes[u]);
    Matrix_3 Q = QuaternionFrame::average(q1, q2).to_rotation_matrix();
    taus[i] = Q.transpose() * (p.at(v) - p.at(u)) * invScale;
    // taus[i] = Eigen::round(taus[i].array());
    i++;
  }
}

/****************************** public ***********************************/

int
Hierarchy::build_hierarchy(const DualGraph& dg,
                           const FrameField& ff,
                           double scale,
                           int depth)
{
  // build base layer;
  layers.clear();
  layers.push_back(HierarchyLayer());
  auto& base_layer = layers[0];
  // vertices
  dg.request_vertices(&(base_layer.vertices));
  int n_vertices = (int)base_layer.vertices.size();
  // quaternion frames
  base_layer.qframes = ff.request_QuaternionFrame();
  // 2-degree graph
  base_layer.edge_list.clear();
  SMatrix<int> am = dg.request_adjacent_matrix();
  Eigen::SparseMatrix<int, Eigen::ColMajor> bm = am;
  am = am * bm;
  for (int i = 0; i < n_vertices; i++) {
    for (SMatrix<int>::InnerIterator jt(am, i); jt; ++jt) {
      double diff = ff[i].diff(ff[jt.col()]);
      base_layer.edge_list.push_back({ i, jt.col() });
    }
  }

  auto& log = DebugTools::Log::get_instance();

  // upcast_matrix.clear();
  upcast_tables.clear();
  FLOGOUT(<< "------------- build hierarchy ------------------")
  while (depth--) {
    if (!build_layer())
      break;
  }
  FLOGOUT(<< "------------- build done ------------------")
  cell_scale = scale;
  return layers.size();
}

void
Hierarchy::optimize()
{
  auto& logger = DebugTools::Log::get_instance();
  logger.flog << "hierarchy optimization scale: " << cell_scale << std::endl;
  int depth = layers.size() - 1;
  int iteration_times = 150;
  optimized_positions = layers[depth].vertices;
  // optimize from top to bottom
  while (depth > 0) {
    logger.tick();
    layers[depth].optimize_positions(
      &optimized_positions, iteration_times, cell_scale);
    positions_down_cast(&optimized_positions, upcast_tables[depth - 1]);
    depth--;
    logger.tick();
    logger.print_time_span_tostd();
  }
  // deal with base layer
  layers[0].optimize_positions(
    &optimized_positions, iteration_times, cell_scale);
}

#ifdef LATTICE_USE_TBB
void
Hierarchy::parallel_optimize()
{
  auto& logger = DebugTools::Log::get_instance();
  logger.flog << "hierarchy optimization scale: " << cell_scale << std::endl;
  size_t depth = layers.size() - 1;
  int iteration_times = 200;
  optimized_positions = layers[depth].vertices;
  // optimize from top to bottom
  while (depth > 0) {
    logger.tick();
    layers[depth].parallel_optimize_positions(
      &optimized_positions, iteration_times, cell_scale);
    positions_down_cast(&optimized_positions, upcast_tables[depth - 1]);
    depth--;
    logger.tick();
    logger.print_time_span_tostd();
    logger.flog << "depth: " << depth + 1
                << " time: " << logger.get_millisec_span() << std::endl;
  }
  // deal with base layer
  logger.tick();
  layers[0].parallel_optimize_positions(
    &optimized_positions, iteration_times, cell_scale);
  logger.tick();
  logger.print_time_span_tostd();
  logger.flog << "depth: " << 0 << " time: " << logger.get_millisec_span()
              << std::endl;
  logger.flog << "optimized vertices: " << optimized_positions.size()
              << std::endl;
}
#endif

void
Hierarchy::dense_optimize(const FrameField* ff)
{
  auto& vertices = dense_layer.vertices;
  vertices = optimized_positions;

  auto& qframes = dense_layer.qframes;
  for (int i = 0; i < vertices.size(); ++i) {
    qframes.push_back(ff->interpolate_qframe(vertices[i]));
  }

  auto& edges = dense_layer.edge_list;
  edges.clear();
  // remove bevel edges
  for (auto& e : layers[0].edge_list) {
    int u = e[0];
    int v = e[1];
    V3d tau =
      cal_tau(vertices[u], vertices[v], qframes[u], qframes[v]) / cell_scale;
    if (is_orthogonal_vector(tau)) {
      edges.push_back(e);
    }
  }

  // insert new vertices
  this->noval_vertices = dense_layer.increase_vertex_density(ff, cell_scale);

  optimized_positions = dense_layer.vertices;
#ifdef LATTICE_USE_TBB
  int iteration_time = 150;
  dense_layer.parallel_optimize_positions(
    &optimized_positions, iteration_time, cell_scale);
#else
  int iteration_time = 70;
  dense_layer.optimize_positions(
    &optimized_positions, iteration_time, cell_scale);
#endif
}
}
// namespace LatticeCore