#include <map>
#include <array>
#include <string>
#include <vector>

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

#include <quadmath.h>

#include "util.h"

// no standard, portable pi until c++20
constexpr double pi = 3.14159265358979323846;

using point_t = std::array<__float128, 2>;

using namespace std;

// p = (d == '+') in original
double angle(point_t &u, point_t &v, point_t &w, bool p=true) {

  auto vu = atan2q(u[1] - v[1], u[0] - v[0]);
  auto vw = atan2q(w[1] - v[1], w[0] - v[0]);
  auto phi = vw - vu;

  if(phi < 0)
    phi += 2.0 * pi;
  if(!p)
    phi = 2.0 * pi - phi;

  return phi;
}

point_t intersect(point_t &a, point_t &b, point_t &c, point_t &d) {
  auto s = (b[0] - a[0]) * (d[1] - c[1]) - (d[0] - c[0]) * (b[1] - a[1]);
  auto x = ((b[0] * a[1] - a[0] * b[1]) * (d[0] - c[0]) - (d[0] * c[1] - c[0] * d[1]) * (b[0] - a[0])) / s;
  auto y = ((b[0] * a[1] - a[0] * b[1]) * (d[1] - c[1]) - (d[0] * c[1] - c[0] * d[1]) * (b[1] - a[1])) / s;
  return {x, y};
}

vector<point_t> optimal_polygon(vector<double> xp, vector<double> yp, double w) {

  assert(xp.size() == yp.size());
  vector<__float128> x(xp.size());
  vector<__float128> y(yp.size());
  for(int k = 0; k < x.size(); k += 1) {
    x[k] = xp[k];
    y[k] = yp[k];
  }

  point_t p_pls = {x[0], y[0] + w};
  point_t l_pls = {x[0], y[0] + w};
  point_t r_pls = {x[1], y[1] + w};

  point_t p_mns = {x[0], y[0] - w};
  point_t l_mns = {x[0], y[0] - w};
  point_t r_mns = {x[1], y[1] - w};

  map<point_t, point_t> s_pls = {{ {x[0], y[0] + w}, {x[1], y[1] + w} }};
  map<point_t, point_t> t_pls = {{ {x[1], y[1] + w}, {x[0], y[0] + w} }};
  map<point_t, point_t> s_mns = {{ {x[0], y[0] - w}, {x[1], y[1] - w} }};
  map<point_t, point_t> t_mns = {{ {x[1], y[1] - w}, {x[0], y[0] - w} }};

  point_t p_i_pls;
  point_t p_i_mns;

  vector<point_t> res;
  uint32_t i = 2;

  while(i < y.size()) {
    // updating ch_plus (convex hull) and ch_minus
    point_t p = {x[i - 1], y[i - 1] + w};
    p_i_pls = {x[i], y[i] + w};
    while (p != p_pls && angle(p_i_pls, p, t_pls.at(p), true) > pi)
      p = t_pls.at(p);
    s_pls[p] = p_i_pls;
    t_pls[p_i_pls] = p;

    p = {x[i - 1], y[i - 1] - w};
    p_i_mns = {x[i], y[i] - w};
    while (p != p_mns && angle(p_i_mns, p, t_mns.at(p), false) > pi)
      p = t_mns.at(p);
    s_mns[p] = p_i_mns;
    t_mns[p_i_mns] = p;

    // check if ch_plus and ch_minus intersect
    if(angle(p_i_pls, l_pls, r_mns, true) < pi) {
      res.push_back(intersect(l_pls, r_mns, p_pls, p_mns));
      p_mns = r_mns;
      point_t t = {x[i - 1], y[i - 1] + w};
      p_pls = intersect(l_pls, r_mns, t, p_i_pls);
      s_pls[p_pls] = p_i_pls;
      t_pls[p_i_pls] = p_pls;
      r_pls = p_i_pls;
      r_mns = p_i_mns;
      l_pls = p_pls;
      l_mns = p_mns;
      while(angle(l_mns, r_pls, s_mns.at(l_mns), false) < pi)
        l_mns = s_mns.at(l_mns);
    } else if (angle(p_i_mns, l_mns, r_pls, false) < pi) {
      res.push_back(intersect(l_mns, r_pls, p_mns, p_pls));
      p_pls = r_pls;
      point_t t = {x[i - 1], y[i - 1] - w};
      p_mns = intersect(l_mns, r_pls, t, p_i_mns);
      s_mns[p_mns] = p_i_mns;
      t_mns[p_i_mns] = p_mns;
      r_mns = p_i_mns;
      r_pls = p_i_pls;
      l_mns = p_mns;
      l_pls = p_pls;
      while(angle(l_pls, r_mns, s_pls.at(l_pls), true) < pi)
        l_pls = s_pls.at(l_pls);
    } else {
      // updating the two separating and supporting lines
      if(angle(p_i_pls, l_mns, r_pls, true) < pi) {
        r_pls = p_i_pls;
        while(angle(p_i_pls, l_mns, s_mns.at(l_mns), true) < pi)
          l_mns = s_mns.at(l_mns);
      }
      if(angle(p_i_mns, l_pls, r_mns, false) < pi) {
        r_mns = p_i_mns;
        while(angle(p_i_mns, l_pls, s_pls.at(l_pls), false) < pi)
          l_pls = s_pls.at(l_pls);
      }
    }
    i += 1;
  }

  // add last change point
  auto a = intersect(l_pls, r_mns, p_pls, p_mns);
  auto b = intersect(l_mns, r_pls, p_mns, p_pls);
  point_t p = {(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0};
  res.push_back(p);

  auto end_a = intersect(p, r_pls, p_i_mns, p_i_pls);
  auto end_b = intersect(p, r_mns, p_i_mns, p_i_pls);
  point_t end = {(end_a[0] + end_b[0]) / 2.0, (end_a[1] + end_b[1]) / 2.0};
  res.push_back(end);

  return res;
}

array<vector<double>, 2> parse(string &file) {

  string line;
  ifstream mstream(file);
  getline(mstream, line);

  array<vector<double>, 2> xy;
  vector<string> cols;
  while(getline(mstream, line)) {
    cols.clear();
    split(line, ',', cols);
    xy[0].push_back(stod(cols[0]));
    xy[1].push_back(stod(cols[1]));
  }

  return xy;
}


int main(int argc, char *argv[]) {

  if(argc < 3) {
    cerr << "usage: plfoptq <input.csv> <w> [precision]\n";
    return 1;
  }

  string p_file(argv[1]);
  string p_w(argv[2]);

  double w = stod(p_w);
  int d = 6;
  if(argc > 3) {
    string p_d(argv[3]);
    d = stoi(p_d);
  }


  cerr << "input: " << p_file << "\n";
  cerr << "    w: " << w << "\n";
  cerr << "    d: " << d << "\n";

  auto input = parse(p_file);
  if(input[0].size() < 3){
    throw runtime_error("input must have at least 3 points");
  }

  auto res = optimal_polygon(input[0], input[1], w);

  auto f = "%1." + to_string(d) + "Qf";
  char buf[128];
  for(const auto &r : res){
    quadmath_snprintf(buf, sizeof(buf), f.c_str(), r[0]);
    cout << buf << ",";
    quadmath_snprintf(buf, sizeof(buf), f.c_str(), r[1]);
    cout << buf << "\n";
  }
  return 0;
}
