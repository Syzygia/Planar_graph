#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <valarray>

#define OUTER_FACE 0

#define dfs_leads \
(int)(edges[edge_ind].leads[0] == prev_vertex)

#define find_leads \
(int)(edges[i].leads[0] == prev_vertex)

#define find_leads_spline \
(int)(edges[i].leads[0] == current_p_ind)

typedef enum Color{white,grey,black} Color;

int BASE_CUR_STEP = 10;

typedef enum Point_type{standard, spline} Point_type;

struct point {
  double x,y;
  std::set<int> faces = {0, 1};
  Color color = white;
  std::vector<int> edges_ind;
  std::vector<int> inserted_edges_ind;
  bool is_in_cycle = false;
  Point_type point_type = standard;

  point() = default;
  point(double x,double y, Point_type pt = standard):x(x), y(y), point_type(pt) {}
  point(const point& tmp):
      x(tmp.x),
      y(tmp.y),
      edges_ind(tmp.edges_ind),
      color(tmp.color),
      point_type(tmp.point_type),
      is_in_cycle(tmp.is_in_cycle),
      inserted_edges_ind(tmp.inserted_edges_ind),
      faces(tmp.faces){}
  point& operator = (point const &tmp){
      x = tmp.x;
      y = tmp.y;
      edges_ind = tmp.edges_ind;
      inserted_edges_ind = tmp.inserted_edges_ind;
      color = tmp.color;
      is_in_cycle = tmp.is_in_cycle;
      point_type = tmp.point_type;
      faces = tmp.faces;
      return  *this;
  }
};


struct face {
  std::unordered_set<int> vertexes;
  std::unordered_set<int> edges;
};

struct edge;
std::vector<face> faces;
std::vector<edge> edges;

struct edge {
  int faces_ind [2] = {0, 1};
  int leads [2];
  bool is_in_cycle = false;
  edge(int first_p, int second_p, int face_a= 0, int face_b = 1) {
      leads[0] = first_p;
      leads[1] = second_p;
      faces_ind[0] = face_a;
      faces_ind[1] = face_b;
      if (face_a != 0 || face_b != 1) {
          faces[face_a].edges.insert(edges.size());
          faces[face_b].edges.insert(edges.size());
      }
  }

};


std::vector<point> points;
std::vector<int> inserted_edges;
std::vector<int> cycle_vertexes;
std::vector<int> cycle_edges;
std::set<int> inserted_points = {};

struct section {
  //using vertex_iter = std::vector<point>::iterator;
  //using edge_iter = std::vector<edge>::iterator;

  std::vector<int> vertexes;
  std::vector<int> edges;
  int com_faces_num = 2;
  std::set<int> com_faces {0, 1};

  //find num of common faces for vertexes
  void get_com_faces_num () {
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
  bool operator < (const section& rhs) const {
      return this->com_faces_num < rhs.com_faces_num;
  }
};
std::vector<section> sections;



//accepts last visited vertex index
bool dfs (int prev_vertex, int prev_edge) {
    cycle_vertexes.push_back(prev_vertex);
    points[prev_vertex].color = grey;
    //cycle over ribs of vertex
    for (auto& edge_ind : points[prev_vertex].edges_ind) {


        //if we checked this edge in undirected graph
        if (edge_ind == prev_edge) {
            continue;
        }

        //find next vertex to call dfs in edge structure

        if (points[edges[edge_ind].leads[dfs_leads]].color) {
            cycle_vertexes.push_back(edges[edge_ind].leads[dfs_leads]);
            return true;
        }
        if (dfs (edges[edge_ind].leads[dfs_leads], edge_ind)) {
            return true;
        }

        points[prev_vertex].color = black;
    }
    return false;
}

int Cycle_length = 0;

//mark edges is_in_cycle and add them to inserted_edges
void get_cycle_edges () {
    auto first = cycle_vertexes.back();
    ++Cycle_length;
    points[first].is_in_cycle = true;
    auto prev = first;
    for (int i = (int) cycle_vertexes.size() - 2; i >= 0; --i) {
        points[cycle_vertexes[i]].is_in_cycle = true;

        //find rib witch leads to prev vertex
        for (auto &j : points[cycle_vertexes[i]].edges_ind) {
            if (edges[j].leads[0] == prev || edges[j].leads[1] == prev) {
                points[cycle_vertexes[i]].inserted_edges_ind.push_back(j);
                points[prev].inserted_edges_ind.push_back(j);
                edges[j].is_in_cycle = true;
                inserted_edges.push_back(j);
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
std::unordered_set <int> segmented;
//contains visited vertexes on a current step of selecting segment
std::unordered_set <int> visited_on_segment;

//function to dfs through segments (excluding cycles)
bool is_path_founded = false;
 void section_dfs (section& elem, const int& edge_ind, int prev_vertex) {
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

void get_sections () {

    int prev_added = -1;
    int ind = 0;
    int contact_vert_num = -1;
    for (auto const &i : edges) {
        contact_vert_num = (points[i.leads[0]].is_in_cycle || points[i.leads[0]].color == black) ?  0 : 1;
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
        section_dfs(elem, ind, i.leads[contact_vert_num == 0? 1 : 0]);
        ++ind;


        sections.push_back(std::move(elem));
    }
}


//inserts segment to  cycle subgraph
using section_iter = std::vector<section>::iterator;


inline int area (point const &a, point const &b,point const &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool intersect_1 (double a, double b, double c, double d) {
    if (a > b)  std::swap (a, b);
    if (c > d)  std::swap (c, d);
    return std::max(a,c) <= std::min(b,d);
}

bool intersect (point const & a, point const & b, point const & c, point const & d) {
    return intersect_1(a.x, b.x, c.x, d.x)
        && intersect_1(a.y, b.y, c.y, d.y)
        && area(a, b, c) * area(a, b, d) <= 0
        && area(c, d, a) * area(c, d, b) <= 0;
}

std::pair<double, double> find_place (int face_ind, double x, double y, double radius, int point_a, int point_b = -1) {
     double x_new = 0;
     double y_new = 0;
     int theta = 0;
     bool is_break = false;
     for (theta = 0; theta <= 360; theta++) {
         x_new = radius * cos(theta);
         y_new = radius * sin(theta);
         is_break = false;
         for (auto const &edge_ind: faces[face_ind].edges) {
            if (intersect(point(x_new, y_new), points[point_a],
                          points[edges[edge_ind].leads[0]], points[edges[edge_ind].leads[1]]) ||
                point_b != -1 && intersect(point(x_new, y_new), points[point_b],
                                           points[edges[edge_ind].leads[0]], points[edges[edge_ind].leads[1]])) {
                is_break = true;
                break;
            }
         }
         if (is_break) {
             continue;
         }
         break;
     }
     if (theta == 361) {
         std::cout << 'can not find place';
         exit;
     }
     return std::make_pair(x_new, y_new);
 }


//connects 2 spline points and gives them edge, and face relation
void add_spline_segment(int prev_p, int new_p) {
    //HUGE ERROR IF 1 FACE WAS CHANGED???
    points[prev_p].edges_ind.push_back(edges.size());
    points[new_p].edges_ind.push_back(edges.size());
    points[new_p].inserted_edges_ind.push_back(edges.size());
    points[prev_p].inserted_edges_ind.push_back(edges.size());

    //next line is kinda pointless?
    //faces[*points[prev_p].faces.rbegin()].vertexes.extract(prev_p);
    points[prev_p].faces.erase(*points[prev_p].faces.rbegin());
    points[prev_p].faces.insert((int) faces.size() - 1);

    auto node = points[new_p].faces.extract(*points[new_p].faces.rbegin());
    node.value() = (int) faces.size() - 1;
    points[new_p].faces.insert(std::move(node));

    faces.back().edges.insert(edges.size());
    faces.back().vertexes.insert(prev_p);
    faces.back().vertexes.insert(new_p);
    faces[0].vertexes.insert(prev_p);
    faces[0].vertexes.insert(new_p);
    inserted_edges.push_back(edges.size());
    edges.emplace_back(prev_p, new_p, 0, faces.size() - 1);
}

 void construct_path(const int& start_p_ind, const int& end_p_ind) {
     //need to check overlap!
     faces.emplace_back();
//     inserted_edges.push_back(edges.size());
//     edges.emplace_back(start_p_ind, points.size(), 0, faces.size() - 1);


     //auto &start_p = points[start_p_ind];
     bool is_on_top = points[start_p_ind].y == 0;

     if (is_on_top) {
        points.emplace_back(points[start_p_ind].x, points[start_p_ind].y + BASE_CUR_STEP, spline);
     }
     else {
         points.emplace_back(points[start_p_ind].x, points[start_p_ind].y - BASE_CUR_STEP, spline);
     }

    points[start_p_ind].faces.insert((int)(faces.size() - 1));

     add_spline_segment(start_p_ind, points.size() -1);
     int prev_vertex = points.size() - 1;
     int direction = start_p_ind % 2 ? -BASE_CUR_STEP : BASE_CUR_STEP;

     //point to run on cycle verge
     auto *current_p = &points[start_p_ind];
     int current_p_ind = start_p_ind;
     point *it_point;
     point *best_next_point;
     double best_x = current_p->x;
     while (true) {

        for (auto const& i : current_p->inserted_edges_ind) {
            it_point = &points[edges[i].leads[find_leads_spline]];

            if (edges[i].leads[find_leads_spline] == end_p_ind) {
                if (it_point->y != current_p->y) {
                    points.emplace_back(best_next_point->x + direction, best_next_point->y, spline);
                    add_spline_segment(prev_vertex, points.size() - 1);
                }
                auto l = *it_point->faces.rbegin();
                add_spline_segment(points.size()- 1, edges[i].leads[find_leads_spline]);
                it_point->faces.insert(l);
                return;
            }
            if ((((it_point->x <= best_x) && (direction == -BASE_CUR_STEP)) ||
                (it_point->x >=  best_x) && (direction == BASE_CUR_STEP)) &&
                edges[i].is_in_cycle) {
                best_x = it_point->x;
                best_next_point = it_point;
                current_p_ind = edges[i].leads[find_leads_spline];
            }
        }
         is_on_top = best_next_point->y == 0;

         if (best_next_point->y != current_p->y){
             points.emplace_back(best_next_point->x + direction, current_p->y, spline);
             add_spline_segment(points.size() - 2, points.size() -1);
             points.emplace_back(best_next_point->x + direction, best_next_point->y, spline);
             add_spline_segment(points.size() - 2, points.size() -1);
             points.emplace_back(best_next_point->x, is_on_top? best_next_point->y + BASE_CUR_STEP :
                         best_next_point->y - BASE_CUR_STEP, spline);
             add_spline_segment(points.size() - 2, points.size() -1);
             //change direction
             direction = direction == BASE_CUR_STEP ?  -BASE_CUR_STEP : BASE_CUR_STEP;
         } else {
             if (is_on_top) {
                 points.emplace_back(best_next_point->x , best_next_point->y + BASE_CUR_STEP, spline);
             }
             else {
                 points.emplace_back(best_next_point->x , best_next_point->y - BASE_CUR_STEP, spline);
             }
             add_spline_segment(prev_vertex, points.size() -1);
         }

        current_p = best_next_point;

        //! mb  incorrect
        faces[*current_p->faces.begin()].vertexes.erase(current_p_ind);
        auto node = current_p->faces.extract(*current_p->faces.begin());
        node.value() = faces.size() - 1;
        current_p->faces.insert(std::move(node));
        faces.back().vertexes.insert(current_p_ind);


        prev_vertex = points.size() - 1;
     }
     //watch out f
 }

 void draw_line(const point &start, const point &end, std::vector<int> const &vertexes, std::vector<int> edges_inds,
                bool is_spline_type, bool is_x_axis, int prev_face) {
    double step_main = 0;
    double step_collateral = 0;
    double direction = 0;
    int count = -1;
    int edge_ind = 0;
    int prev_vertex_ind = -1;
    for (auto const &j: vertexes) {
        if (prev_vertex_ind != -1) {
            points[j].inserted_edges_ind.push_back(edges_inds[edge_ind - 1]);
        }
        prev_vertex_ind = j;
        if (edge_ind < edges_inds.size()) {
            points[j].inserted_edges_ind.push_back(edges_inds[edge_ind++]);
        }
        ++count;
        if (points[j].is_in_cycle || *vertexes.begin() == j || vertexes.back() == j) {
            points[j].faces.insert({prev_face, (int)faces.size()});
            continue;
        } else {
            points[j].faces = {prev_face, (int)faces.size()};
        if (!is_x_axis) {
            step_main = (end.x - start.x) / (double) vertexes.size();
            step_collateral = (end.y - start.y) / (double) vertexes.size();
            direction = start.y != 0? BASE_CUR_STEP / (double) vertexes.size() :
                -BASE_CUR_STEP / (double) vertexes.size();
            points[j].x = start.x + step_main * count;
            points[j].y = start.y + (is_spline_type ? direction : 0) + step_collateral * count;
        } else {
            step_main = (end.y - start.y) / (double) vertexes.size();
            step_collateral = (end.x - start.x) / (double) vertexes.size();
            direction = start.x == 0? BASE_CUR_STEP / (double) vertexes.size():
                -BASE_CUR_STEP / (double) vertexes.size();
            points[j].y = start.y + step_main * count;
            points[j].x = start.x + (is_spline_type ? direction : 0) + step_collateral * count;
            }
        }
    }
}



void insert_segment (section_iter &section_iter) {
    auto face_to_sep = std::max_element(section_iter->com_faces.begin(), section_iter->com_faces.end());
    if (*face_to_sep == OUTER_FACE) {
        std::copy(section_iter->vertexes.begin(), section_iter->vertexes.end(),
                  std::inserter(inserted_points, inserted_points.end()));
        return construct_path(*section_iter->vertexes.begin(), *(++section_iter->vertexes.begin()));
    }
    auto &prev_face = faces[*face_to_sep];
    face new_face;
    double k_coeff = (points[*(section_iter->vertexes.begin())].y - points[(section_iter->vertexes.back())].y) /
        (points[*(section_iter->vertexes.begin())].x - points[(section_iter->vertexes.back())].x);
    double b_coeff = points[*(section_iter->vertexes.begin())].y -
        k_coeff * points[*(section_iter->vertexes.begin())].x;
    bool is_horisontal = points[*(section_iter->vertexes.begin())].x == points[section_iter->vertexes.back()].x;
    double x_border = points[*(section_iter->vertexes.begin())].x;
    int face_num = 0;


    //buffer for iteration + erase
    auto prev_vertexes = prev_face.vertexes;

    bool is_contains = false;
    //move a part of vertexes to new face and assign new face relation in them
    for (auto const &point_ind: prev_vertexes) {
        //due to double type nature statement bellow may not function as intended
        if ((is_horisontal && points[point_ind].x == x_border) ||
            (points[point_ind].y == k_coeff * points[point_ind].x + b_coeff)) {
            points[point_ind].faces.insert((int) faces.size());
            auto prev_vertex = point_ind;
            for (auto const &edge_ind: points[point_ind].inserted_edges_ind) {
                if ((is_horisontal && points[edges[edge_ind].leads[dfs_leads]].x <= x_border) ||
                        (points[edges[edge_ind].leads[dfs_leads]].y >=
                        k_coeff * points[edges[edge_ind].leads[dfs_leads]].x + b_coeff)) {
                    face_num = edges[edge_ind].faces_ind[0] == *face_to_sep ? 0 : 1;
                    if (edges[edge_ind].faces_ind[face_num] != *face_to_sep) {
                        continue;
                    }
                    prev_face.edges.erase(edge_ind);
                    new_face.edges.insert(edge_ind);
                    edges[edge_ind].faces_ind[face_num] = (int) faces.size();
                }
            }
        }
        if ((is_horisontal && points[point_ind].x >= x_border) ||
        (points[point_ind].y <= k_coeff * points[point_ind].x + b_coeff)) {
            continue;
        }


        new_face.vertexes.insert(point_ind);
        auto node = points[point_ind].faces.extract(*face_to_sep);
        node.value() = (int) faces.size();
        points[point_ind].faces.insert(std::move(node));
        prev_face.vertexes.erase(point_ind);
        //remove all edges connected to point witch is no longer in this face
        for (auto const &edge_ind: points[point_ind].inserted_edges_ind) {
            face_num = edges[edge_ind].faces_ind[0] == *face_to_sep ? 0 : 1;
            if (edges[edge_ind].faces_ind[face_num] != *face_to_sep) {
                continue;
            }
            edges[edge_ind].faces_ind[face_num] = (int) faces.size();
            prev_face.edges.erase(edge_ind);
            new_face.edges.insert(edge_ind);
        }


    }

    prev_face.edges.insert(section_iter->edges.begin(), section_iter->edges.end());
    inserted_edges.insert(inserted_edges.end(), section_iter->edges.begin(), section_iter->edges.end());
    new_face.vertexes.insert(section_iter->vertexes.begin(), section_iter->vertexes.end());
    prev_face.vertexes.insert(section_iter->vertexes.begin(), section_iter->vertexes.end());
    for (auto const &i: section_iter->edges) {
        //mb wrong!!
        edges[i].faces_ind[0] = *face_to_sep;
        edges[i].faces_ind[1] = (int)faces.size();
    }
    if (section_iter->vertexes.size() > 2) {
        //need to check whether contact points always on edge of vector!!!
        auto contact_a = points[*section_iter->vertexes.begin()];
        int current_p_ind = *section_iter->vertexes.begin();
        auto contact_b = points[section_iter->vertexes.back()];
        double step = (contact_b.x - contact_a.x) / (double) section_iter->vertexes.size();
        bool is_spline_type = false;
        if (contact_a.x == contact_b.x) {
            for (auto const &i : contact_a.inserted_edges_ind) {
                if ((edges[i].is_in_cycle || inserted_points.contains(edges[i].leads[find_leads_spline]))  &&
                points[edges[i].leads[find_leads_spline]].x == contact_a.x) {
                    is_spline_type = true;
                    break;
                }
            }
            draw_line(contact_a, contact_b, section_iter->vertexes, section_iter->edges, is_spline_type,
                      true, *face_to_sep);
        } else {
            for (auto const &i : contact_a.inserted_edges_ind) {
                if ((edges[i].is_in_cycle || inserted_points.contains(edges[i].leads[find_leads_spline])) &&
                points[edges[i].leads[find_leads_spline]].y == contact_a.y) {
                    is_spline_type = true;
                    break;
                }
            }
            draw_line(contact_a, contact_b, section_iter->vertexes, section_iter->edges, is_spline_type,
                      false, *face_to_sep);
        }
    }


    //?
    new_face.edges.insert(section_iter->edges.begin(), section_iter->edges.end());
    std::copy(section_iter->vertexes.begin(), section_iter->vertexes.end(),
              std::inserter(inserted_points, inserted_points.end()));
    faces.push_back(new_face);

 }

 void read_data() {
    int size;
    std::fstream input;
    input.open("tests/1.txt", std::ios_base::in);
    input >> size;
    points.resize(size);

    int a_point, b_point;
    while(input >> a_point >> b_point) {
        points[a_point].edges_ind.push_back(edges.size());
        points[b_point].edges_ind.push_back(edges.size());
        edge ed (a_point, b_point);
        edges.push_back(ed);
     }
}


//
void init_faces() {
    face outer_face;
    face inner_face;
    bool first = false;
    for (auto const &i: cycle_vertexes) {
        if (points[i].is_in_cycle &&
        first)  {
            outer_face.vertexes.insert(i);
            inner_face.vertexes.insert(i);
        }
        first = true;
    }

    for (auto const &i: cycle_edges) {
            outer_face.edges.insert(i);
            inner_face.edges.insert(i);
    }
    faces.push_back(std::move(outer_face));
    faces.push_back(std::move(inner_face));
}


void assign_pos_cycle () {
    int j = 0;
    double last_x = -BASE_CUR_STEP;
    double last_y = 0;
    for (int i = (int)cycle_vertexes.size() -1; i > 0; --i) {
        if (i != cycle_vertexes.size() -1 && cycle_vertexes[i] == cycle_vertexes[cycle_vertexes.size() -1]) {
            break;
        }

        // draw a rectangle with cycle vertexes

        if (j < Cycle_length  / 2 ) {
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

void make_dot_file () {
    std::ofstream dot_file ("1.dot");
    dot_file << "graph G {" << std::endl;
    int j = 0;
    for (auto const &i: points) {
        if (i.point_type != spline) {
            dot_file << j++ << " [pos = \"" << i.x << "," <<i.y << "!\"]" << std::endl;
        } else {
            dot_file << j++ << " [pos = \"" << i.x << "," <<i.y << "!\"" << " style=\"invis\"" << "]" << std::endl;
        }
    }
    for (auto const &i: inserted_edges) {
        dot_file << edges[i].leads[0] << "--" << edges[i].leads[1] << std::endl;
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
            if (!faces[j].vertexes.contains(vertex_ind)) {
                std::cout << "face correctness failed in vertex " << vertex_ind << " and face " << j << std::endl;
            }
        }
        ++vertex_ind;
    }
    for (auto const &i: inserted_edges) {
        for (auto const &j: edges[i].faces_ind) {
            if (!faces[j].edges.contains(i)) {
                std::cout << "face correctness failed in edge " << i << " and face " << j << std::endl;
            }
        }
    }

}

int main() {
    read_data();
    if (!dfs(0, -1)) {
        std::cout << "no cycle";
        exit(1);
    }
    get_cycle_edges();
    assign_pos_cycle();
    init_faces();
    get_sections();

    while (!sections.empty()) {

        //find segment with lowest num of intersected faces
        auto b_sect_iter = std::min_element(sections.begin(), sections.end());
        if (b_sect_iter == sections.end()) {
            std::cout << "poop";
        }
        std::cout << b_sect_iter->com_faces_num << std::endl;
        if (b_sect_iter->com_faces_num == 0) {
            std::cout << "segment with 0 intersected faces found, non-planar graph";
            exit(1);
        }

        insert_segment(b_sect_iter);

        sections.erase(b_sect_iter);
        //recalculate faces values of segments after last insertion
        for(auto &i : sections) {
            i.get_com_faces_num();
        }
    }
    make_dot_file();

    check_face_correctness();
    return 0;
}
