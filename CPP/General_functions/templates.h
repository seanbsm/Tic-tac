#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <vector>
#include <numeric>

/* Copied from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes */
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), (size_t) 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  //std::stable_sort(idx.begin(), idx.end(),
  /* The use of this algorithm is a guaranteed case of unequal entries,
   * -> use the faster std:sort() */
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

#endif // TEMPLATES_H