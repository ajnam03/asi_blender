#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/mesh_segmentation.h>
#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;


std::string createSegmentMeshFilename(std::size_t iteration, std::size_t segmentIndex)
{
    std::stringstream ss;
    ss << "iteration_" << iteration << "_segment_" << segmentIndex << ".off";
    return ss.str();
}


int main(int argc, char* filepath[])

    //  // Load your mesh here ---> MAYBE ADD THESE LINES
    // Mesh mesh;

    // // Triangulate the faces of the mesh
    // CGAL::Polygon_mesh_processing::triangulate_faces(mesh);


{   
    // read in mesh
    const std::string filename = filepath[1];
    Mesh mesh;
    if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
    std::cerr << "Invalid input." << std::endl;
    return 1;
    }
  
//   // remesh AN_EDIT
//   double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.04;
//   unsigned int nb_iter = 3;
//   std::cout << "Split border...";
//   std::vector<edge_descriptor> border;
//   PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
//   PMP::split_long_edges(border, target_edge_length, mesh);
//   std::cout << "done." << std::endl;
//   std::cout << "Start remeshing of " << filename
//     << " (" << num_faces(mesh) << " faces)..." << std::endl;
//   PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
//                            CGAL::parameters::number_of_iterations(nb_iter)
//                                             .protect_constraints(true)); //i.e. protect border, here
//   std::cout << "Remeshing done." << std::endl;
//AN EDIT --
  // SM mesh; 
  // Define ranges of parameter options to iterate through
  // std::vector<double> sdfValueParameters = {0.1, 0.2, 0.3};
  // std::vector<double> segmentationOptions = {0.1, 0.2, 0.3};
  std::vector<double> sdfConeAngles = {60};
  std::vector<double> sdfRayNumbers = {25};
  std::vector<double> segNumClusters = {2,3};
  std::vector<double> segSmoothLambdas = {0.05,0.1, 0.5, 0.55};
//   std::vector<double> segNumClusters = {2, 3};
//   std::vector<double> segSmoothLambdas = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55};


  for (const double sdfConeAng : sdfConeAngles)
  {
    for (const double sdfRayNum : sdfRayNumbers)
    {
        for (const double segNumClust : segNumClusters)
        {
            for (const double segSmoothL : segSmoothLambdas)
            {
                std::cout << "cluster #: " << segNumClust << std::endl; //print iteration details
                std::cout << "lambda smooth level: " << segSmoothL << std::endl;

                // smooth mesh
                const unsigned int nb_iterations = 5;
                const double time = 0.0001;
                std::set<Mesh::Vertex_index> constrained_vertices;
                for(Mesh::Vertex_index v : vertices(mesh))
                {
                if(is_border(v, mesh))
                    constrained_vertices.insert(v);
                }
                std::cout << "Constraining: " << constrained_vertices.size() << " border vertices" << std::endl;
                CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);
                std::cout << "Smoothing shape... (" << nb_iterations << " iterations)" << std::endl;
                PMP::smooth_shape(mesh, time, CGAL::parameters::number_of_iterations(nb_iterations)
                                                                .vertex_is_constrained_map(vcmap));


                typedef Mesh::Property_map<face_descriptor,double> Facet_double_map;
                Facet_double_map sdf_property_map;
                sdf_property_map = mesh.add_property_map<face_descriptor,double>("f:sdf").first;
                  CGAL::sdf_values(mesh, sdf_property_map, sdfConeAng/180 * CGAL_PI, sdfRayNum, true);
  
                // create a property-map for segment-ids
                typedef Mesh::Property_map<face_descriptor, std::size_t> Facet_int_map;
                Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor,std::size_t>("f:sid").first;;
                // segment the mesh using default parameters for number of levels, and smoothing lambda
                // Any other scalar values can be used instead of using SDF values computed using the CGAL function

                std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, segNumClust, segSmoothL);
                typedef CGAL::Face_filtered_graph<Mesh> Filtered_graph;
                //print area of each segment and then put it in a Mesh and print it in an OFF file if number of segments > 0 and <4
                
                std::cout << "num of segments : " << number_of_segments << std::endl; //print number of segments

                if (number_of_segments == 2)
                {
                    Filtered_graph segment_mesh(mesh);
                    for(std::size_t id = 0; id < number_of_segments; ++id)
                    {
                    segment_mesh.set_selected_faces(id, segment_property_map);
                    std::cout << "Segment "<<id<<"'s area is : "<<CGAL::Polygon_mesh_processing::area(segment_mesh)<<std::endl;
                    Mesh out;
                    CGAL::copy_face_graph(segment_mesh, out);
                    std::ostringstream oss;
                    std::string s = filepath[1];
                    std::string delimiter = ".";
                    std::string new_name = s.substr(0, s.find(delimiter)); // 
                    oss << new_name << "_seg" << id << ".ply";
                    CGAL::IO::write_polygon_mesh(oss.str(), out, CGAL::parameters::stream_precision(17));
                    }
                    goto end;
                }
             
           

            }
        }
    }
  }
end:
std::cout << "segmentation routine complete" << std::endl; 
}
