#pragma once
#include "eigenDefs.h"
// #define use_nanoflann
#ifdef use_nanoflann
#include "nanoflann.hpp"
class KD3d
{
  // typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<> Kdtree_t;
public:
  KD3d(std::vector<V3d>* points) {}
};
#else
#include "gishi523.hpp"
struct KDPoint3d
{
  V3d point;
  static const int DIM = 3;
  KDPoint3d()
    : point(0, 0, 0)
  {}
  KDPoint3d(const V3d& p)
    : point(p)
  {}
  double operator[](int id) const { return point(id); }
};

class KD3d : public kdt::KDTree<KDPoint3d>
{
public:
  KD3d(const ColMatrix<double>& point)
    : KDTree<KDPoint3d>()
  {
    std::vector<KDPoint3d> kdp;
    int n = point.cols();
    for (int i = 0; i < n; ++i) {
      kdp.push_back(KDPoint3d(V3d(point.col(i))));
    }
    build(kdp);
  }

  KD3d(const std::vector<V3d>& point)
  {
    std::vector<KDPoint3d> kdp;
    int n = point.size();
    for (int i = 0; i < n; ++i) {
      kdp.push_back(KDPoint3d(point[i]));
    }
    build(kdp);
  }

  const V3d& get_point(int id) const { return points_[id].point; }
  
  std::vector<int> radiusSearch(const V3d& query, double radius) const
  {
    return KDTree<KDPoint3d>::radiusSearch(KDPoint3d(query), radius);
  }

  std::vector<int> kknnSearch(const V3d& query, int k) const
  {
	  return KDTree<KDPoint3d>::knnSearch(KDPoint3d(query), k);
  }
  int nearestSearch(const V3d& query, double* minDist = nullptr) const
  {
	  return KDTree<KDPoint3d>::nnSearch(KDPoint3d(query), minDist);
  }

};
#endif