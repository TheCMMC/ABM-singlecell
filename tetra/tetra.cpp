#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Epick::Point_3 MyPoint;
typedef CGAL::Triangulation_vertex_base_with_info_3<MyPoint, K> VV;
typedef CGAL::Triangulation_data_structure_3<VV> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
using namespace std;

extern "C" {
  void cfd_setup(char *);
  void cfd_query(double, double, double, double *, double *, double *);
};

class Tetra {
  Triangulation T;
  bool initialized = false;
public:
  bool filled(){return initialized;}
  void fill_triangulation(string &fname) {
    // read in lines of six double precision floating point numbers:
    // x y z u v w
    // representing a point (x,y,z) and its associated velocity (u, v, w).
    std::ifstream data_file(fname.c_str());
    std::vector<std::pair<MyPoint,Point>> points;
    std::string line;
    std::cerr<<"Reading "<<fname<<"... ";
    while (getline(data_file, line)) {
      if (line[0] == '%' || line[0] == '\n' || line[0] == '\0') continue;
      std::istringstream iss(line);
      double x, y, z, u, v, w;
      iss>>x;iss>>y;iss>>z;iss>>u;iss>>v;iss>>w;
      points.push_back(std::make_pair(Point(x,y,z),MyPoint(u,v,w)));
    }
    std::cerr<<" read "<<points.size()<<" points.\n";
    T.insert(points.begin(), points.end());
    initialized = true;
  }

  void estimate_velocity(double x, double y, double z,
			 double &u, double &v, double &w) {
    assert(T.is_valid());
    Locate_type lt;
    int li, lj;
    Point p(x,y,z);
    clock_t t = clock();
    Cell_handle c = T.locate(p, lt, li, lj);
    //    t = clock()-t; std::cerr<<"(Locate in "<<double(t)/CLOCKS_PER_SEC<<" secs): ";

    //  If query lies inside the affine hull of the points, the k-face
    //(finite or infinite) that contains query in its interior is
    //returned, by means of the cell returned together with lt, which is
    //set to the locate type of the query (VERTEX, EDGE, FACET, CELL, or
    //OUTSIDE_CONVEX_HULL if the cell is infinite and query lies
    //strictly in it) and two indices li and lj that specify the k-face
    //of the cell containing query.  If the k-face is a cell, li and lj
    //have no meaning; if it is a facet (resp. vertex), li gives the
    //index of the facet (resp. vertex) and lj has no meaning; if it is
    //and edge, li and lj give the indices of its vertices.  If the
    //point query lies outside the affine hull of the points, which can
    //happen in case of degenerate dimensions, lt is set to
    //OUTSIDE_AFFINE_HULL, and the cell returned has no meaning. As a
    //particular case, if there is no finite vertex yet in the
    //triangulation, lt is set to OUTSIDE_AFFINE_HULL and locate returns
    //the default constructed handle.
    //
    // Need clear up how to access vertices by index to fill in the
    // switch table below, though most random points will be in CELLs
    // and properly handled
    
    MyPoint vel;
    u = v = w = 0;
    switch (lt) {
    case Triangulation::VERTEX:
      //      std::cerr<<"("<<x<<", "<<y<<", "<<z<<") is a Vertex\n";
      /*
	vel = c->vertex(li)->info();
	u = vel.x(); v = vel.y(); w = vel.z();
      */
      break;
    case Triangulation::EDGE:
      //      std::cerr<<"("<<x<<", "<<y<<", "<<z<<") is in an Edge\n";
      /*
	for (int i = 0; i < 2; i++) {
        vel = c->edges(li)->vertex(i)->info();
        u += vel.x();
        v += vel.y();
        w += vel.z();
	}
      */
      break;
    case Triangulation::FACET:
      //      std::cerr<<"("<<x<<", "<<y<<", "<<z<<") is on a Facet\n";
      /*
	for (int i = 0; i < 3; i++) {
        vel = c->select_facet(li)->vertex(i)->info();
        u += vel.x();
        v += vel.y();
        w += vel.z();
	}
	u/=3; v/=3; w/=3;
      */
      break;
    case Triangulation::CELL:
      // Average the vertex velocities as a very rough approximation.
      // Better would be to interpolate properly.
      for (int i = 0; i < 4; i++) {
        vel = c->vertex(i)->info();
        u += vel.x();
        v += vel.y();
        w += vel.z();
      }
      u/=4; v/=4; w/=4;
      break;
    case Triangulation::OUTSIDE_CONVEX_HULL:
      //      std::cerr<<"("<<x<<", "<<y<<", "<<z<<") is outside the Convex Hull\n";
      break;
    case Triangulation::OUTSIDE_AFFINE_HULL:
      //      std::cerr<<"("<<x<<", "<<y<<", "<<z<<") is outside the Affine Hull\n";
      break;
    default:
      std::cerr<<"ERROR: "<<__FILE__<<":"<<__LINE__<<" in estimate_velocity\n";
    }
  }
};

static Tetra T;
void cfd_setup (char * fname) {
  string f = std::string(fname);
  if (!T.filled()) T.fill_triangulation(f);
}
void cfd_query (double x, double y, double z, double * u, double * v, double * w) {
  T.estimate_velocity(x, y, z, *u, *v, *w);
}
