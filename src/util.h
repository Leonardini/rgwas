#include <string>
#include <vector>
#include <cstdint>
#include <sstream>

using namespace std;

inline vector<string> split(const string &str, char delim) {
  vector<string> elems{};
  stringstream ss;
  ss.str(str);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

inline void split(const string &str, char delim, vector<string> &elems) {
  elems.clear();
  stringstream ss;
  ss.str(str);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

// https://stackoverflow.com/questions/15301885/calculate-value-of-n-choose-k
inline uint64_t choose(uint64_t n, uint64_t k) {
  if(k == 0)
    return 1;
  unsigned long long c;
  if(__builtin_umulll_overflow(n, choose(n - 1, k - 1), &c))
    throw overflow_error("nCk overflow");
  return c / k;
}

template<typename T>
inline uint32_t popcnt(T prof, size_t blocks) {
  size_t size = sizeof(prof[0]) * blocks / sizeof(uint64_t);
  auto *data = (uint64_t *) prof;
  uint32_t popcnt = 0;
  for(int k = 0; k < size; k += 1)
    popcnt += __builtin_popcountl(data[k]);
  return popcnt;
}

inline uint32_t popcnt(const vector<uint64_t> &prof) {
  uint32_t popcnt = 0;
  for(const auto k : prof)
    popcnt += __builtin_popcountl(k);
  return popcnt;
}
