#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <numeric>

#define OUTER_FACE 0

#define dfs_leads \
(int)(edges[edge_ind].leads[0] == prev_vertex)

#define find_leads \
(int)(edges[i].leads[0] == prev_vertex)

typedef enum Color{white,grey,black} Color;

int BASE_CUR_STEP = 1;

struct point {
  double x,y;
  std::unordered_set<int> faces = {0, 1};
  Color color = white;
  std::vector<int> edges_ind;
  bool is_in_cycle = false;
  point() = default;
  point(double x,double y):x(x), y(y) {}
  point(const point& tmp):
      x(tmp.x),
      y(tmp.y),
      edges_ind(tmp.edges_ind),
      color(tmp.color),
      is_in_cycle(tmp.is_in_cycle){}
  point& operator = (point const &tmp){
      x = tmp.x;
      y = tmp.y;
      edges_ind = tmp.edges_ind;
      color = tmp.color;
      is_in_cycle = tmp.is_in_cycle;
      return *this;
  }
};

struct edge {
  int faces [2] = {0, 1};
  int leads [2];
  bool is_in_cycle = false;
  edge(int first_p, int second_p) {
      leads[0] = first_p;
      leads[1] = second_p;
  }
};


struct face {
  std::unordered_set<int> vertexes;
  std::unordered_set<int> edges;
};

std::vector<face> faces;
std::vector<point> points;
std::vector<edge> edges;
std::vector<int> cycle_vertexes;
std::vector<int> cycle_edges;

struct section {
  //using vertex_iter = std::vector<point>::iterator;
  //using edge_iter = std::vector<edge>::iterator;

  std::unordered_set<int> vertexes;
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
          if (!points[i].is_in_cycle) {
              continue;
          }
          for (auto const &j : points[i].faces) {
              buf2.insert(j);
          }
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
std::unordered_set <int> segmented;
//contains visited vertexes on a current step of selecting segment
std::unordered_set <int> visited_on_segment;

//function to dfs through segments (excluding cycles)
 void section_dfs (section& elem, const int& edge_ind, int prev_vertex) {
    if (visited_on_segment.count(prev_vertex)) {
        return;
    }
    elem.vertexes.insert(prev_vertex);
    visited_on_segment.insert(prev_vertex);
    if (points[prev_vertex].is_in_cycle && elem.vertexes.size() >= 2) {
        return;
    }
    for (auto &i : points[prev_vertex].edges_ind) {
        if (!edges[i].is_in_cycle && !segmented.count(i)) {
            segmented.insert(i);
            elem.edges.push_back(i);
            section_dfs(elem, i, edges[i].leads[find_leads]);
        }
    }

}

void get_sections () {

    int prev_added = -1;
    int ind = 0;
    for (auto const &i : edges) {
        if (i.is_in_cycle || segmented.count(ind)) {
            ++ind;
            continue;
        }
        section elem;
        elem.vertexes.insert(i.leads[0]);
        elem.edges.push_back(ind);
        segmented.insert(ind);
        visited_on_segment.clear();
        section_dfs(elem, ind, i.leads[1]);
        ++ind;


        sections.push_back(std::move(elem));
    }
}


//inserts segment to  cycle subgraph
using section_iter = std::vector<section>::iterator;

 void construct_path(const int& start_p_ind, const int& end_p_ind) {
     //need to check overlap!
     edges.emplace_back(start_p_ind, points.size());
     auto const& start_p = points[start_p_ind];
     bool is_on_top = start_p.y == 0;

     if (is_on_top) {
        points.emplace_back(start_p.x, start_p.y + BASE_CUR_STEP);
     }
     else {
         points.emplace_back(start_p.x, start_p.y - BASE_CUR_STEP);
     }
     auto &last_p = points.back();
     last_p.edges_ind.push_back(edges.size());
     int prev_vertex = points.size() - 1;
     int direction = start_p_ind % 2 ? -BASE_CUR_STEP : BASE_CUR_STEP;
     face new_face;

     //point to run on cycle verge
     auto &current_p = start_p;
     auto  &it_point = const_cast<point &>(start_p);
     while (true) {
        for (auto const& i : current_p.edges_ind) {
            it_point = points[edges[i].leads[find_leads]];
            if ((it_point.x == current_p.x ||
                it_point.x ==  current_p.x + direction) &&
                it_point.is_in_cycle) {
                is_on_top = it_point.y == 0;
                if (is_on_top) {
                    points.emplace_back(it_point.x , it_point.y + BASE_CUR_STEP);
                }
                else {
                    points.emplace_back(it_point.x , it_point.y - BASE_CUR_STEP);
                }
                if (it_point.y != current_p.y){
                   points.emplace_back(it_point.x + direction, it_point.y);

                   //change direction
                   direction = direction == BASE_CUR_STEP ?  -BASE_CUR_STEP : BASE_CUR_STEP;
                }
                current_p = points[points.size()];
            }
        }
        break;
     }
 }

 void one_edge_face_separation (section_iter section_iter) {
     auto face_to_sep = std::max_element(section_iter->com_faces.begin(), section_iter->com_faces.end());
     if (*face_to_sep == OUTER_FACE) {
        construct_path(*section_iter->vertexes.begin(), *(++section_iter->vertexes.begin()));
     }
     auto &prev_face = faces[*face_to_sep];
     face new_face;
     int ind = 0;
     bool is_equal = false;
     double x_border = -100;
     x_border = points[*(section_iter->vertexes.begin())].x;
//     if (x_border == -100) {
//         std::cout << "no contact vertexes in segment";
//         exit(-1);
//     }
//    if (points[*(iter->vertexes.begin())].x == points[*(++iter->vertexes.begin())].x) {
//        is_equal = true;
//    }

    //buffer for iteration + erase
    auto prev_vertexes = prev_face.vertexes;
    //move a part of vertexes to new face and assign new face relation in them
    for (auto const &point_ind: prev_vertexes) {

         if (points[point_ind].x >= x_border && !section_iter->vertexes.contains(point_ind)) {
             continue;
         }

         new_face.vertexes.insert(point_ind);

         if (!section_iter->vertexes.contains(point_ind)) {
             auto node = points[point_ind].faces.extract(*face_to_sep);
             node.value() = (int) faces.size();
             points[point_ind].faces.insert(std::move(node));
             prev_face.vertexes.erase(point_ind);
             //remove all edges connected to point witch is no longer in this face
             for (auto const &edge_ind: points[point_ind].edges_ind) {
                 // subject to change!!!!
                 edges[edge_ind].faces[1] = faces.size();
                 //double erase??
                 prev_face.edges.erase(edge_ind);
             }

         }
         else {
             points[point_ind].faces.insert((int) faces.size());
         }
    }


    prev_face.edges.insert(section_iter->edges.begin(), section_iter->edges.end());
    edges[section_iter->edges[0]].faces[0] = *face_to_sep;
    edges[section_iter->edges[0]].faces[1] = faces.size();
    //!!!!!!!!!!! potential memory leak
    new_face.edges = std::move(std::unordered_set(section_iter->edges.begin(), section_iter->edges.end()));

    //change face relation in ribs of inserted segment
    for (auto const &i: new_face.vertexes) {
         for (auto const &j : points[i].edges_ind) {
             //!
             if (points[edges[j].leads[0]].faces.contains(faces.size()) &&
                 points[edges[j].leads[1]].faces.contains(faces.size())) {
                 // this is kinda wrong!!!!!!!!!!!!!
                 edges[j].faces[1] = faces.size();
                 new_face.edges.insert(j);
        //                 if (!section_iter->vertexes.contains(i)) {
        //                    prev_face.edges.erase(j);
        //                 }
             }
        }
    }
    faces.push_back(new_face);
 }


void insert_segment (section_iter iter) {

    if (iter->vertexes.size() == 2) {
        one_edge_face_separation(iter);
    }
    // insert segmented part in graph, in new face
    //TODO OUTER FACE, insert of complex segments



 }

 void read_data() {
    int size;
    std::fstream input;
    input.open("3.txt", std::ios_base::in);
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
        if (points[i].is_in_cycle) {
            outer_face.edges.insert(i);
            inner_face.edges.insert(i);
        }
    }
    faces.push_back(std::move(outer_face));
    faces.push_back(std::move(inner_face));
}


void assign_pos_cycle () {
    int j = 0;
    double last_x = -1;
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
        dot_file << j++ << " [pos = \"" << i.x << "," <<i.y << "!\"]" << std::endl;
    }
    for (auto const &i: edges) {
        dot_file << i.leads[0] << "--" << i.leads[1] << std::endl;
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
            if (edges[j].faces[0] != face_ind && edges[j].faces[1] != face_ind) {
                std::cout << "face correctness failed in face " << face_ind << " and edge " << j << std::endl;
            }
        }
        ++face_ind;
    }

    int vertex_ind = 0;
    int edge_ind = 0;
    for (auto const &i: points) {
        for (auto const &j: i.faces) {
            if (!faces[j].vertexes.contains(vertex_ind)) {
                std::cout << "face correctness failed in vertex " << vertex_ind << " and face " << j << std::endl;
            }
        }
        ++vertex_ind;
    }
    for (auto const &i: edges) {
        for (auto const &j: i.faces) {
            if (!faces[j].edges.contains(edge_ind)) {
                std::cout << "face correctness failed in edge " << edge_ind << " and face " << j << std::endl;
            }
        }
        ++edge_ind;
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
    get_sections();
    init_faces();

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
