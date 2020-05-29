#pragma once
#include <stack>
#include <vector>

class MergeFindSet
{
protected:
  mutable std::vector<int> marks;

public:
  MergeFindSet(int size)
  {
    marks.resize(size);
    for (int i = 0; i < size; ++i) {
      marks[i] = i;
    }
  }

  size_t size() const { return marks.size(); }

  int find(int id) const
  {
    if (marks[id] == id)
      return id;
    std::stack<int> s;
    while (!s.empty())
      s.pop();
    int k = id;
    while (marks[k] != k) {
      s.push(k);
      k = marks[k];
    }
    int res = k;
    while (!s.empty()) {
      k = s.top();
      s.pop();
      marks[k] = res;
    }
    return res;
  }

  void merge(int x, int y)
  {
    if (find(y) != find(x))
      marks[y] = marks[x];
  }

  int operator[](int id) const { return find(id); }
};