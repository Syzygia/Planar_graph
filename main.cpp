#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <valarray>
#include <cmath>
#include <random>

#define OUTER_FACE 0

#define dfs_leads \
(int)(edges[edge_ind].leads[0] == prev_vertex)

#define find_leads \
(int)(edges[i].leads[0] == prev_vertex)

#define find_leads_spline \
(int)(edges[i].leads[0] == current_p_ind)

typedef enum Color { white, grey, black } Color;

int BASE_CUR_STEP = 100;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//template <class T>
class vertexes_in_face: public std::vector<int> {
     public:
        int &operator [](int i) {
            i = i < 0 ? i * (-1) - 1 : i;
            i = i >= this->size() ? i - (int)this->size() : i;
            return std::vector<int>::operator[](i);
        }
};

const char *__asan_default_options() {
    return "strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1";
}

typedef enum Point_type { standard, hidden } Point_type;

struct point {
  double x, y;
  std::set<int> faces = {0, 1};
  Color color = white;
  std::vector<int> edges_ind;
  std::vector<int> drawn_edges_ind;
  bool is_in_cycle = false;
  Point_type point_type = standard;

  point() = default;
  point(double x, double y, Point_type pt = standard) : x(x), y(y), point_type(pt){}
  point(const point &tmp) = default;
  point &operator=(point const &tmp) = default;
};

std::vector<point> points;

struct face {
  vertexes_in_face vertexes;
  std::unordered_set<int> edges;
  std::pair<double, double> get_center_mass() const {
      double sum_x = 0, sum_y = 0;
      for (auto const &i: vertexes) {
          sum_x += points[i].x;
          sum_y += points[i].y;
      }
      return std::make_pair<double, double>(sum_x / (double) vertexes.size(), sum_y / (double) vertexes.size());
  }
  face &operator = (face const &a) = default;

  void erase_vertex(int const &vertex_num) {
      auto place = std::find(vertexes.begin(), vertexes.end(), vertex_num);
      vertexes.erase(place);
  }
};

std::vector<face> faces;

struct edge {
  int faces_ind[2] = {0, 1};
  int leads[2];
  bool is_in_cycle = false;
  //bool is_inserted = false;
  edge(int first_p, int second_p) {
      leads[0] = first_p;
      leads[1] = second_p;
  }

};


class Edges {
 public:
    std::vector<edge> &get() {
        return edges;
    }

    std::set<int> &get_inserted_ed() {
        return inserted_edges;
    }

    edge &operator [](int i) {
        return edges[i];
    }
    //inserts edge between point a and b
    void insert_edge(int p_a, int p_b, int f_a = 0, int f_b = 1, int ed_ind = -1) {
     if (ed_ind == - 1) {
         ed_ind = (int) edges.size();
         auto &it = edges.emplace_back(p_a, p_b);
     }

     edges[ed_ind].faces_ind[0] = f_a;
     edges[ed_ind].faces_ind[1] = f_b;

     points[p_a].drawn_edges_ind.push_back(ed_ind);
     points[p_b].drawn_edges_ind.push_back(ed_ind);
     if (f_a != 0 || f_b != 1) {
         faces[f_a].edges.insert(ed_ind);
         faces[f_b].edges.insert(ed_ind);
     }
     inserted_edges.insert(ed_ind);
    }

 private:
    std::vector<edge> edges;
    std::set<int> inserted_edges;
};

Edges edges;

std::vector<int> cycle_vertexes;
std::vector<int> cycle_edges;
std::set<int> inserted_points = {};

//holds segments (path) witch are going to be inserted in graph
struct section {
  //using vertex_iter = std::vector<point>::iterator;
  //using edge_iter = std::vector<edge>::iterator;

  std::vector<int> vertexes;
  std::vector<int> edges;
  int com_faces_num = 2;
  std::set<int> com_faces{0, 1};

  //find num of common faces for vertexes
  void get_com_faces_num() {
      std::set<int> buf1;
      for (int i = 0; i < faces.size(); ++i) {
          buf1.insert(buf1.end(), i);
      }
      std::set<int> buf2;
      std::set<int> buf3;
      for (auto const &i : vertexes) {
          buf2.clear();
          //!
          if (!points[i].is_in_cycle && !inserted_points.contains(i)) {
              continue;
          }
          buf2.insert(points[i].faces.begin(), points[i].faces.end());
          std::set_intersection(buf1.begin(), buf1.end(), buf2.begin(), buf2.end(),
                                std::inserter(buf3, buf3.begin()));
          buf1 = buf3;
          buf3.clear();
      }
      com_faces_num = buf1.size();
      com_faces = buf1;
  }
  bool operator<(const section &rhs) const {
      return this->com_faces_num < rhs.com_faces_num;
  }
};
std::vector<section> sections;

//accepts last visited vertex index
bool dfs(int prev_vertex, int prev_edge) {
    cycle_vertexes.push_back(prev_vertex);
    points[prev_vertex].color = grey;
    //cycle over ribs of vertex
    for (auto &edge_ind : points[prev_vertex].edges_ind) {


        //if we checked this edge in undirected graph
        if (edge_ind == prev_edge) {
            continue;
        }

        //find next vertex to call dfs in edge structure

        if (points[edges[edge_ind].leads[dfs_leads]].color) {
            cycle_vertexes.push_back(edges[edge_ind].leads[dfs_leads]);
            return true;
        }
        if (dfs(edges[edge_ind].leads[dfs_leads], edge_ind)) {
            return true;
        }

        points[prev_vertex].color = black;
    }
    return false;
}

int Cycle_length = 0;

//mark edges is_in_cycle and add them to inserted_edges
void get_cycle_edges() {
    auto first = cycle_vertexes.back();
    ++Cycle_length;
    points[first].is_in_cycle = true;
    auto prev = first;
    for (int i = (int) cycle_vertexes.size() - 2; i >= 0; --i) {
        points[cycle_vertexes[i]].is_in_cycle = true;

        //find rib witch leads to prev vertex
        for (auto &j : points[cycle_vertexes[i]].edges_ind) {
            if (edges[j].leads[0] == prev || edges[j].leads[1] == prev) {
                edges.insert_edge(cycle_vertexes[i], prev, 0, 1, j);
                edges[j].is_in_cycle = true;
                cycle_edges.push_back(j);
                break;
            }
        }
        if (cycle_vertexes[i] == first) {
            break;
        }
        ++Cycle_length;
        prev = cycle_vertexes[i];
    }

}

//contains all visited edges when finding segments (prevents copy of segments)
std::unordered_set<int> segmented;
//contains visited vertexes on a current step of selecting segment
std::unordered_set<int> visited_on_segment;

//function to dfs through segments (excluding cycles)
bool is_path_founded = false;
void section_dfs(section &elem, const int &edge_ind, int prev_vertex) {
    if (visited_on_segment.contains(prev_vertex) || is_path_founded) {
        return;
    }
    elem.vertexes.push_back(prev_vertex);
    visited_on_segment.insert(prev_vertex);
    if ((points[prev_vertex].is_in_cycle || points[prev_vertex].color == black) &&
        elem.vertexes.size() >= 2) {
        is_path_founded = true;
        return;
    }
    points[prev_vertex].color = black;
    for (auto &i : points[prev_vertex].edges_ind) {
        if (!edges[i].is_in_cycle && !segmented.count(i) && !is_path_founded) {
            segmented.insert(i);
            elem.edges.push_back(i);
            section_dfs(elem, i, edges[i].leads[find_leads]);
        }
    }

}

void get_sections() {

    int prev_added = -1;
    int ind = 0;
    int contact_vert_num = -1;
    for (auto const &i : edges.get()) {
        contact_vert_num = (points[i.leads[0]].is_in_cycle || points[i.leads[0]].color == black) ? 0 : 1;
        if (i.is_in_cycle || segmented.count(ind) ||
            !(points[i.leads[contact_vert_num]].is_in_cycle ||
                points[i.leads[contact_vert_num]].color == black)) {
            ++ind;
            continue;
        }
        section elem;
        elem.vertexes.push_back(i.leads[contact_vert_num]);
        elem.edges.push_back(ind);
        segmented.insert(ind);
        visited_on_segment.clear();
        is_path_founded = false;
        section_dfs(elem, ind, i.leads[contact_vert_num == 0 ? 1 : 0]);
        ++ind;

        sections.push_back(std::move(elem));
    }
}


void embed_anchor_point_new_face(int const &face_ind, int const &point_ind) {
    faces[face_ind].vertexes.push_back(point_ind);
    points[point_ind].faces.insert(face_ind);
}

//manages int-between points faces and their edges during face separation
void embed_interim_point_new_face(int const &face_ind, int const &prev_face_ind, int const &point_ind) {
    faces[face_ind].vertexes.push_back(point_ind);
    faces[prev_face_ind].erase_vertex(point_ind);
    auto node = points[point_ind].faces.extract(prev_face_ind);
    node.value() = face_ind;
    points[point_ind].faces.insert(std::move(node));

    int face_num;
    //remove all edges connected to point witch is no longer in this face
    for (auto const &edge_ind: points[point_ind].drawn_edges_ind) {
        face_num = edges[edge_ind].faces_ind[0] == prev_face_ind ? 0 : 1;
        if (edges[edge_ind].faces_ind[face_num] != prev_face_ind) {
            continue;
        }
        edges[edge_ind].faces_ind[face_num] = (int) faces.size() - 1;
        faces[prev_face_ind].edges.erase(edge_ind);
        faces[face_ind].edges.insert(edge_ind);
    }
}

//inserts segment to  cycle subgraph
using section_iter = std::vector<section>::iterator;

const double EPS = 1E-9;

inline double det(double a, double b, double c, double d) {
    return a * d - b * c;
}

inline bool between(double a, double b, double c) {
    return std::min(a, b) <= c + EPS && c <= std::max(a, b) + EPS;
}

inline bool intersect_1(double a, double b, double c, double d) {
    if (a > b) std::swap(a, b);
    if (c > d) std::swap(c, d);
    return std::max(a, c) < std::min(b, d);
}

bool a_equal(double a, double b)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}
bool intersect(point const &a, point const &b, point const &c, point const &d) {
    double A1 = a.y - b.y, B1 = b.x - a.x, C1 = -A1 * a.x - B1 * a.y;
    double A2 = c.y - d.y, B2 = d.x - c.x, C2 = -A2 * c.x - B2 * c.y;
    double zn = det(A1, B1, A2, B2);

    if (!a_equal(zn, 0)) {
        double x = -det(C1, B1, C2, B2) * 1. / zn;
        double y = -det(A1, C1, A2, C2) * 1. / zn;
//        if ((x == c.x && y == c.y) || (x == d.x && y == d.y)) {
//            return false;
//        }
//      pls dont be ugly
        if ((a_equal(x,b.x) && a_equal(y, b.y)) || (a_equal(x, a.x) && a_equal(y, a.y))) {
            return false;
        }
        return between(a.x, b.x, x) && between(a.y, b.y, y)
            && between(c.x, d.x, x) && between(c.y, d.y, y);
    } else
        return det(A1, C1, A2, C2) == 0 && det(B1, C1, B2, C2) == 0
            && intersect_1(a.x, b.x, c.x, d.x)
            && intersect_1(a.y, b.y, c.y, d.y);
}


std::tuple<bool, double, double> check_edge_placement(int const &face_ind, double const &theta,
                                                      double const &x, double const &y, double radius,
                                                      int point_a ,int point_b = -1) {
    double x_new = x + radius * cos((double)theta / 180 * M_PI);
    double y_new = y + radius * sin((double)theta / 180 * M_PI);
    bool is_break = false;
    point temp_p;
    temp_p.x = x_new;
    temp_p.y = y_new;
    for (auto const &edge_ind: faces[face_ind].edges) {
        if (!edges.get_inserted_ed().contains(edge_ind)) {
            continue;
        }
        if (point_b != -1 && intersect(temp_p, points[point_b],
                                       points[edges[edge_ind].leads[0]], points[edges[edge_ind].leads[1]]) ||
            intersect(temp_p, points[point_a],
                      points[edges[edge_ind].leads[0]], points[edges[edge_ind].leads[1]])) {
            is_break = true;
            break;
        }
    }
    return std::make_tuple(is_break, x_new, y_new);
}

//tries to place a point on a given face with 360 sector and given radius
std::pair<double, double> find_place(int face_ind, double x, double y, double radius, int point_a , int point_b = -1,
                                     double start_angle = 361, double end_angle = 361) {
    double x_new = 0;
    double y_new = 0;
    int theta = 0;
    bool is_break = false;
    radius = fabs(radius);
//    double temp = std::max(start_angle, end_angle);
//    start_angle = std::min(start_angle, end_angle);
//    end_angle = temp;
    if (start_angle == 361 || end_angle == 361) {
        start_angle = 0;
        end_angle = 360;
    } else {
        std::tie(is_break, x_new, y_new) = check_edge_placement(face_ind,
                                                                ((int)(ceil(start_angle + (start_angle + end_angle) / 2)) % 360),
                                                                x, y, radius, point_a, point_b);
        if (!is_break) {
            return std::make_pair(x_new, y_new);
        }
    }
    bool is_place_found = false;
    for (theta = ceil(start_angle); theta != ceil(end_angle); theta = (theta + 1) % 360) {
        is_break = false;
        std::tie(is_break, x_new, y_new) =
            check_edge_placement(face_ind, theta, x, y, radius, point_a, point_b);
        if (is_break) {
            continue;
        }

        // is this always work?
        is_place_found = true;
        break;
    }
    if (!is_place_found) {
        std::cout << "can not find place";
        std::exit(1);
    }
    return std::make_pair(x_new, y_new);
}

//attempts to place (many) points in a path between anchor points
void draw_line(const point &start, const point &end, std::vector<int> const &vertexes, std::vector<int> const &edges_inds,
               bool is_spline_type, int prev_face) {
    double step_main = 0;
    double step_collateral = 0;
    int count = -1;
    int edge_ind = 0;
    int prev_vertex_ind = -1;
    double radius = 0;
    for (auto const &j: vertexes) {
        prev_vertex_ind = j;
        ++count;
        if (points[j].is_in_cycle || *vertexes.begin() == j || vertexes.back() == j) {
            points[j].faces.insert({prev_face, (int) faces.size() - 1});
            //inserted points?
            continue;
        } else {
            points[j].faces = {prev_face, (int) faces.size() - 1};
            step_main = (end.x - start.x) / (double) vertexes.size();
            step_collateral = (end.y - start.y) / (double) vertexes.size();
            std::pair<double, double> coord;
            if (is_spline_type) {
                if (fabs(step_collateral) < fabs(step_main)) {
                    radius = fabs(step_collateral) == 0 ?  fabs(step_main) : fabs(step_collateral);
                } else {
                    radius = fabs(step_main) == 0? fabs(step_collateral): fabs(step_main);
                }
//                std::cout << "step_coratel " << fabs(step_collateral) << " step_main " << fabs(step_main) << " radius " <<
//                radius << std::endl;
                if (j != vertexes[vertexes.size() - 2]) {
                    coord = find_place(prev_face, start.x + step_main * count, start.y + step_collateral * count,
                                       radius, vertexes[count - 1]);
                } else {

                    coord = find_place(prev_face, start.x + step_main * count, start.y + step_collateral * count,
                                       radius,vertexes[count - 1], vertexes[count + 1]);
                }
                points[j].x = coord.first;
                points[j].y = coord.second;
            } else {
                points[j].x = start.x + step_main * count;
                points[j].y = start.y + step_collateral * count;
            }
        }
    }
}

double get_vec_angle(point const &a, point const &b, bool is_outer_face) {
//    double dot, det;
//    dot = a.x * b.x + a.y * b.y;
//    det = a.x * b.y - a.y * b.x;
//    return atan2(det, dot) * 180 / M_PI;
    double dot, det;
    //double vec_b_x = b.x - a.x;
    double vec_b_x = 1;
    double vec_b_y = 0;
    double vec_d_x = b.x - a.x;
    double vec_d_y = b.y - a.y;

    dot = vec_b_x * vec_d_x + vec_b_y * vec_d_y;
    det = vec_b_x * vec_d_y - vec_b_y * vec_d_x;
    return atan2(det, dot) * 180 / M_PI;
}

double get_distance(point &a, point &b) {
    double x, y;
    x = (b.x - a.x) / 2;
    y = (b.y - a.y) / 2;
    return x <= y ? x : y;
}

//tries to place point near another point in the given face
void attach_point_point (int const &fc_ind, point &p_to_place, int const &prev_vertex,
                         int const & left_vertex, int const &right_vertex, int const &destination) {
    std::pair<double, double> coord;

    //this is bad! cause of incorrect angle and works dor only given test
    coord = find_place(fc_ind, points[prev_vertex].x, points[prev_vertex].y,
                       fc_ind == OUTER_FACE ? 1 : get_distance(points[prev_vertex], points[destination]),
                       prev_vertex, -1,
                       get_vec_angle(points[prev_vertex],points[left_vertex], fc_ind == OUTER_FACE),
                       get_vec_angle(points[prev_vertex],points[right_vertex], fc_ind == OUTER_FACE));
    p_to_place.x = coord.first;
    p_to_place.y = coord.second;
    //warning! if face is really small - wont work!!

}


double polygonArea(std::vector<double> const &X, std::vector<double> const &Y)
{
    // Initialize area
    double area = 0.0;
    int n = X.size();
    int j = n - 1;
    for (int i = 0; i < n; i++)
    {
        area += (X[j] + X[i]) * (Y[j] - Y[i]);
        j = i;  // j is previous vertex to i
    }

    return (area / 2.0);
}

std::vector<int> get_path_to_anchor_point (int const &fc_ind, int const &start_p, int const & end_p) {
    auto const &fc = faces[fc_ind];
    auto it = fc.vertexes.begin();
    bool is_start_found = false;
    std::vector<int> path;
   while (true) {
       if (*it == start_p && !is_start_found) {
            is_start_found = true;
       }
       if (is_start_found) {
           path.push_back(*it);
       }
       if (*it == end_p) {
           return path;
       }
       ++it;
       if (it == fc.vertexes.end()) {
           it = fc.vertexes.begin();
       }
   }
}

// attempts to draw an edge in not straight line fashion
void draw_along(int const &prev_face, std::vector<int> const &path,vertexes_in_face &vt_in_fc) {
    auto &new_face = *faces.rbegin();
    int prev_point = *path.begin();

    int i = (int)(std::find(vt_in_fc.begin(), vt_in_fc.end(), prev_point) - vt_in_fc.begin());
    for (auto const &point_ind: path) {
        if (point_ind == *path.begin() || point_ind == path.back()) {
            embed_anchor_point_new_face((int)faces.size() -1, point_ind);
            if (prev_face != OUTER_FACE || point_ind == path.back()) {
                ++i;
                continue;
            }
        } else {
            embed_interim_point_new_face((int)faces.size() -1, prev_face,point_ind);
        }
            //!!!!!!!!!!!!!!!!!!!!!!
        point &new_p = points.emplace_back(0, 0, hidden);
        points.back().faces = {prev_face, (int) faces.size()};
        new_p.faces = {prev_face, (int) faces.size()};
        if (prev_face != OUTER_FACE) {
            attach_point_point(prev_face, new_p, point_ind, vt_in_fc[i + 1], vt_in_fc[i - 1],
                               path.back());
        } else {
            attach_point_point(prev_face, new_p, point_ind, vt_in_fc[i - 1], vt_in_fc[i + 1],
                               path.back());
        }
        edges.insert_edge(prev_point, (int) points.size() - 1, prev_face, (int) faces.size() - 1);
        inserted_points.insert((int) points.size() - 1);
        if (point_ind == *(++path.rbegin())) {
            edges.insert_edge((int) points.size() - 1, path.back(), prev_face, (int) faces.size() - 1);
        }
        prev_point = (int) points.size() - 1;
        ++i;
    }
}


void anti_clock_wise_face_insertion (std::vector<int> const &vertexes, int face_ind) {
    //care for position of inserted elem, it might be the end
    bool is_clock_wise_segment = false;
    bool have_found = false;
    for (auto i = faces[face_ind].vertexes.begin(); i <faces[face_ind].vertexes.end(); i++) {
        if (!have_found) {
            if (*i == *vertexes.begin()) {
                have_found = true;
                is_clock_wise_segment = false;
            }
            if (*i == *vertexes.rbegin()) {
                have_found = true;
                is_clock_wise_segment = true;
            }
            if (have_found) {
                if (*(i + 1) == (is_clock_wise_segment? *vertexes.begin() : *vertexes.rbegin())) {
                    if (is_clock_wise_segment) {
                        faces[face_ind].vertexes.insert(++i, vertexes.rbegin()  +1, vertexes.rend()  - 1);
                    }
                    else {
                        faces[face_ind].vertexes.insert(++i, vertexes.begin() +1, vertexes.end() - 1);
                    }
                    break;
                }
            }
        }
        if (have_found && *i == (is_clock_wise_segment? *vertexes.begin(): *vertexes.rbegin())) {
            if (is_clock_wise_segment) {
                faces[face_ind].vertexes.insert(++i, vertexes.begin() +1, vertexes.end() - 1);
            }
            else {
                faces[face_ind].vertexes.insert(++i, vertexes.rbegin()  +1, vertexes.rend()  - 1);
            }
            break;
        }

    }
}


bool is_horisontal_angle_of_rotation(std::vector<int> const &vertexes) {
    std::vector<double> x_es(vertexes.size());
    std::vector<double> y_es(vertexes.size());
    int j = 0;
    for (auto i: vertexes) {
        x_es[j] = points[i].x;
        y_es[j++] = points[i].y;
    }
    auto area = polygonArea(x_es, y_es);
    return area > 0;
}

//inserts complex segment in faces preserving anti-clockwise rotation
void change_face_relation_for_path_segment (std::vector<int> const &vertexes, int prev_face, int new_face) {
    anti_clock_wise_face_insertion(vertexes, prev_face);
    //anti_clock_wise_face_insertion(vertexes, new_face);
    if (is_horisontal_angle_of_rotation(vertexes)) {
        //faces[new_face].vertexes = {vertexes.rbegin(), vertexes.rend()};
        faces[new_face].vertexes.clear();
        faces[new_face].vertexes.insert(faces[new_face].vertexes.begin(), vertexes.rbegin(), vertexes.rend());
    } else {
        faces[new_face].vertexes.clear();
        faces[new_face].vertexes.insert(faces[new_face].vertexes.begin(), vertexes.begin(), vertexes.end());
//        faces[new_face].vertexes = {vertexes.begin(), vertexes.end()};
    }
}


void insert_segment(section_iter &section_iter) {
    auto face_to_sep = *std::max_element(section_iter->com_faces.begin(), section_iter->com_faces.end());
    if (face_to_sep == OUTER_FACE) {
//        std::copy(section_iter->vertexes.begin(), section_iter->vertexes.end(),
//                  std::inserter(inserted_points, inserted_points.end()));
        //return outer_face_pathing(*section_iter->vertexes.begin(), *(++section_iter->vertexes.begin()));
    }
    auto &new_face = faces.emplace_back();
    if (section_iter->vertexes.size() <= 2) {
        std::vector<int> path_to_draw;
        for (auto const &edge_ind :faces[face_to_sep].edges) {
            if (face_to_sep == OUTER_FACE || intersect(points[edges[*section_iter->edges.begin()].leads[0]],
                          points[edges[*section_iter->edges.begin()].leads[1]],
                          points[edges[edge_ind].leads[0]],
                          points[edges[edge_ind].leads[1]])) {
                path_to_draw = get_path_to_anchor_point(face_to_sep, *section_iter->vertexes.begin(),
                                                        *section_iter->vertexes.rbegin());
                draw_along(face_to_sep,path_to_draw, faces[face_to_sep].vertexes);
                return;
            }
        }
    }
    if (section_iter->vertexes.size() > 2) {
        //need to check whether contact points always on edge of vector!!!
        auto contact_a = points[*section_iter->vertexes.begin()];
        int current_p_ind = *section_iter->vertexes.begin();
        auto contact_b = points[section_iter->vertexes.back()];
        double step = (contact_b.x - contact_a.x) / (double) section_iter->vertexes.size();
        bool is_spline_type = false;

        for (auto const &i : contact_a.drawn_edges_ind) {
            if (edges[i].leads[find_leads_spline] == section_iter->vertexes.back()) {
                is_spline_type = true;
                break;
            }
        }
        draw_line(contact_a, contact_b, section_iter->vertexes, section_iter->edges, is_spline_type, face_to_sep);

        //
        //force test !
        if (section_iter->vertexes.size() == 3) {
            points[section_iter->vertexes[1]].y = -1.5;
        }
        //
    }

    double k_coeff = (points[*(section_iter->vertexes.begin())].y - points[(section_iter->vertexes.back())].y) /
        (points[*(section_iter->vertexes.begin())].x - points[(section_iter->vertexes.back())].x);
    double b_coeff = points[*(section_iter->vertexes.begin())].y -
        k_coeff * points[*(section_iter->vertexes.begin())].x;
    bool is_horisontal = points[*(section_iter->vertexes.begin())].x == points[section_iter->vertexes.back()].x;
    double x_border = points[*(section_iter->vertexes.begin())].x;
    int face_num = 0;


    //buffer for iteration + erase
    auto prev_vertexes = faces[face_to_sep].vertexes;
    bool is_anchor_adjacent = false;

    //determine whether or not anchor points are neighbors
    auto prev_vertex = section_iter->vertexes[0];
    for (auto edge_ind : points[prev_vertex].drawn_edges_ind) {
        if (edges[edge_ind].leads[dfs_leads] == *section_iter->vertexes.rbegin()) {
            is_anchor_adjacent = true;
            break;
        }
    }

    //move a part of vertexes to new face and assign new face relation in them
    for (auto const &point_ind: prev_vertexes) {
        //due to double type nature statement bellow may not function as intended
        //if anchor point basically
        if ((is_horisontal && points[point_ind].x == x_border) ||
            (!is_horisontal && points[point_ind].y == k_coeff * points[point_ind].x + b_coeff)) {
            embed_anchor_point_new_face((int)faces.size() - 1, point_ind);
            prev_vertex = point_ind;
            for (auto const &edge_ind: points[point_ind].drawn_edges_ind) {
                //moves ribs witch are no longer belong to old face
                if ((is_anchor_adjacent &&
                (is_horisontal && points[edges[edge_ind].leads[dfs_leads]].x == x_border) ||
                    (!is_horisontal && points[edges[edge_ind].leads[dfs_leads]].y ==
                        k_coeff * points[edges[edge_ind].leads[dfs_leads]].x + b_coeff))
                || (!is_anchor_adjacent &&
                ((is_horisontal && points[edges[edge_ind].leads[dfs_leads]].x <= x_border) ||
                    (!is_horisontal && points[edges[edge_ind].leads[dfs_leads]].y >=
                        k_coeff * points[edges[edge_ind].leads[dfs_leads]].x + b_coeff)))) {
                    face_num = edges[edge_ind].faces_ind[0] == face_to_sep ? 0 : 1;
                    if (edges[edge_ind].faces_ind[face_num] != face_to_sep) {
                        continue;
                    }
                    faces[face_to_sep].edges.erase(edge_ind);
                    new_face.edges.insert(edge_ind);
                    edges[edge_ind].faces_ind[face_num] = (int) faces.size() - 1;
                }
            }
        }
        if ((is_horisontal && points[point_ind].x >= x_border) ||
            (!is_horisontal && points[point_ind].y <= k_coeff * points[point_ind].x + b_coeff)) {
            continue;
        }
        // special case when points on the left should not be moved
        if (is_anchor_adjacent) {
            continue;
        }

        embed_interim_point_new_face((int) faces.size() - 1, face_to_sep, point_ind);

    }
    if (section_iter->vertexes.size() > 2) {
        change_face_relation_for_path_segment(section_iter->vertexes, face_to_sep, (int) faces.size() -1);
    }

//    new_face.vertexes.push_back(section_iter->vertexes.begin(), section_iter->vertexes.end());
//    faces[face_to_sep].vertexes.push_back(section_iter->vertexes.begin(), section_iter->vertexes.end());
    for (auto const &i: section_iter->edges) {
        edges.insert_edge(edges[i].leads[0], edges[i].leads[1], face_to_sep, (int) faces.size() - 1, i);
    }

    //?
    std::copy(section_iter->vertexes.begin(), section_iter->vertexes.end(),
              std::inserter(inserted_points, inserted_points.end()));
}

void read_data() {
    int size;
    std::fstream input;
    input.open("tests/100nodes.txt", std::ios_base::in);
    input >> size;
    points.resize(size);

    int a_point, b_point;
    while (input >> a_point >> b_point) {
        points[a_point].edges_ind.push_back((int)edges.get().size());
        points[b_point].edges_ind.push_back((int)edges.get().size());
        edges.get().emplace_back(a_point, b_point);
    }
}

//
void init_faces(std::vector<face> &fc, std::vector<int> const &cyc_vet) {
    face outer_face;
    face inner_face;
    bool first = false;
    for (auto const &i: cyc_vet) {
        if (points[i].is_in_cycle &&
            first) {
            outer_face.vertexes.push_back(i);
            inner_face.vertexes.push_back(i);
        }
        first = true;
    }

    for (auto const &i: cycle_edges) {
        outer_face.edges.insert(i);
        inner_face.edges.insert(i);
    }
    fc.push_back(std::move(outer_face));
    fc.push_back(std::move(inner_face));
}

void assign_pos_cycle() {
    int j = 0;
    double last_x = -BASE_CUR_STEP;
    double last_y = 0;
    for (int i = (int) cycle_vertexes.size() - 1; i > 0; --i) {
        if (i != cycle_vertexes.size() - 1 && cycle_vertexes[i] == cycle_vertexes[cycle_vertexes.size() - 1]) {
            break;
        }

        // draw a rectangle with cycle vertexes

        if (j < Cycle_length / 2) {
            points[cycle_vertexes[i]].x = last_x + BASE_CUR_STEP;
            points[cycle_vertexes[i]].y = last_y;
            last_x = points[cycle_vertexes[i]].x;
        } else {
            //!
            if (points[cycle_vertexes[i + 1]].y == 0) {
                last_x += BASE_CUR_STEP;
            }
            points[cycle_vertexes[i]].x = last_x - BASE_CUR_STEP;
            points[cycle_vertexes[i]].y = last_y - BASE_CUR_STEP;
            last_x = points[cycle_vertexes[i]].x;
        }
        ++j;
    }
}

void make_dot_file() {
    std::ofstream dot_file("1.dot");
    dot_file << "graph G {" << std::endl;
    int j = 0;
    for (auto const &i: points) {
        if (i.point_type != hidden) {
            dot_file << j++ << " [pos = \"" << i.x << "," << i.y << "!\"" << " width=0.2, height=0.2 " << "]" << std::endl;
        } else {
            dot_file << j++ << " [pos = \"" << i.x << "," << i.y << "!\"" << " width=0, height=0, margin=0, style=\"invis\"" << "]" << std::endl;
        }
    }
    for (auto const &i: edges.get_inserted_ed()) {
        if (points[edges[i].leads[0]].point_type == hidden && points[edges[i].leads[1]].point_type == hidden) {
            dot_file << edges[i].leads[0] << "--" << edges[i].leads[1] << " [headclip=\"false\"]"<< std::endl;
        } else {
            dot_file << edges[i].leads[0] << "--" << edges[i].leads[1] << std::endl;
        }
    }
    dot_file << "}";
    dot_file.close();
}

void check_face_correctness() {
    int face_ind = 0;
    for (auto const &i: faces) {
        for (auto const &j: i.vertexes) {
            if (!points[j].faces.contains(face_ind)) {
                std::cout << "face correctness failed in face " << face_ind << " and vertex " << j << std::endl;
            }
        }

        for (auto const &j: i.edges) {
            if (edges[j].faces_ind[0] != face_ind && edges[j].faces_ind[1] != face_ind) {
                std::cout << "face correctness failed in face " << face_ind << " and edge " << j << std::endl;
            }
        }
        ++face_ind;
    }

    int vertex_ind = 0;
    for (auto const &i: points) {
        for (auto const &j: i.faces) {
//            if (!faces[j].vertexes.contains(vertex_ind)) {
//                std::cout << "face correctness failed in vertex " << vertex_ind << " and face " << j << std::endl;
//            }
        }
        ++vertex_ind;
    }
    for (auto const &i: edges.get_inserted_ed()) {
        for (auto const &j: edges[i].faces_ind) {
            if (!faces[j].edges.contains(i)) {
                std::cout << "face correctness failed in edge " << i << " and face " << j << std::endl;
            }
        }
    }

}

template <typename I>
I random_element(I begin, I end, std::mt19937 &gen)
{
    const unsigned long n = std::distance(begin, end);
    //const unsigned long divisor = (RAND_MAX + 1) / n;

    //unsigned long k;
    //do { k = std::rand() / divisor; } while (k >= n);
    std::uniform_int_distribution dst(1, (int)n - 1);

    std::advance(begin, dst(gen));
    return begin;
}

int point_placed_num = 4;
void place_line_test(std::fstream &input, int between_a, int between_b, int amount) {
    int times = point_placed_num + amount;
    bool is_first = true;
    for (point_placed_num; point_placed_num < times; point_placed_num++) {
        if (is_first) {
            input << between_a << " " << point_placed_num << std::endl;
            is_first = false;
        } else {
            input << point_placed_num - 1 << " " << point_placed_num << std::endl;
            if (point_placed_num == times -1) {
                input << point_placed_num << " " << between_b << std::endl;
            }
        }
    }
}

void generate_tests(int n) {
    std::fstream input;
    input.open("tests/100nodes.txt", std::ios_base::out);
    input << 100 << std::endl;
    input << "0 1\n"
             "1 2\n"
             "2 3\n"
             "3 0\n";

    std::vector<int> cycle_contact_points = {2, 1, 3};
    //int left = 0, right = 2;
    int prev_point = 0;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd());

    std::vector<int> line(20);
    std::iota (line.begin(), line.end(), point_placed_num);
    place_line_test(input, 0, 2, 20);
    auto it = random_element(line.begin(), line.end(), gen);

    std::vector<int> line_1(20);
    std::iota (line_1.begin(), line_1.end(), point_placed_num);
    place_line_test(input, *it, 1, 20);

    std::uniform_int_distribution l_or_r(0, 1);
    int left, right = 0;
    int is_right = 0;
    int num_of_lines = 0;
    unsigned long distance = 0;
    while (num_of_lines < 5) {
        is_right =  l_or_r(gen);
        if (is_right) {
            auto buf_it = it;
            it = random_element(buf_it, line.end(), gen);
        } else {
            it = random_element(line.begin(), it, gen);
        }
        auto it_1 = random_element(line_1.begin(), line_1.end(), gen);
        distance = std::distance(line_1.begin(), it_1);
        left = *it;
        right = *it_1;
        line = std::move(line_1);
        //line = line_1;
        line_1 = std::vector<int> (20);
        std::iota (line_1.begin(), line_1.end(), point_placed_num);
        place_line_test(input, left, right, 20);
        it = line.begin();
        std::advance(it, distance);
        ++num_of_lines;
    }
    exit(0);
    
}

int main() {
    //generate_tests(0);
    read_data();
    if (!dfs(0, -1)) {
        std::cout << "no cycle";
        exit(1);
    }
    get_cycle_edges();
    assign_pos_cycle();
    init_faces(faces, cycle_vertexes);
    get_sections();

    auto is_dirst = false;
    section lowest;
    while (!sections.empty()) {

        //find segment with lowest num of intersected faces
        //auto b_sect_iter = std::min_element(sections.begin(), sections.end());
        section_iter b_sect_iter;
        lowest.com_faces_num = 1000;
        for (auto it=sections.begin(); it != sections.end(); it ++) {
            if (*it < lowest) {
                lowest = *it;
                b_sect_iter = it;
            }
        }
        if (b_sect_iter == sections.end()) {
            std::cout << "poop";
        }
//        //force test
        if (b_sect_iter->vertexes.back() == 5 && !is_dirst) {
            is_dirst = true;
            b_sect_iter = sections.begin() + ((int)sections.size() - 3);

        }
        std::cout << b_sect_iter->com_faces_num << std::endl;
        if (b_sect_iter->com_faces_num == 0) {
            std::cout << "segment with 0 intersected faces found, non-planar graph";
            exit(1);
        }

        insert_segment(b_sect_iter);

        sections.erase(b_sect_iter);
        //recalculate faces values of segments after last insertion
        for (auto &i : sections) {
            i.get_com_faces_num();
        }
    }
    make_dot_file();

    check_face_correctness();
    return 0;
}
