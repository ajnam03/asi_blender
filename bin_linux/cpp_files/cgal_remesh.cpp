#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <iostream>
#include <string>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};
int main(int argc, char* filepath[])
{   
    // read in mesh
    const std::string filename = filepath[1];
    Mesh mesh;
    if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
    std::cerr << "Invalid input." << std::endl;
    return 1;
    }
    // remesh
    double target_edge_length = 0.04;
    unsigned int nb_iter = 3;
    std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
    std::cout << "done." << std::endl;
    std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;
    PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
                            CGAL::parameters::number_of_iterations(nb_iter)
                                            .protect_constraints(true)); //i.e. protect border, here
    // save mesh
    std::ostringstream oss;
    std::string s = filepath[1];
    std::string delimiter = ".";
    std::string new_name = s.substr(0, s.find(delimiter)); // t
    oss << new_name << "_remesh.ply";
    CGAL::IO::write_polygon_mesh(oss.str(), mesh, CGAL::parameters::stream_precision(17));
    std::cout << "Remeshing done." << std::endl;
    return 0;
}
