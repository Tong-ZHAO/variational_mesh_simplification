#ifndef TYPES_H
#define TYPES_H

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

// std
#include <vector>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point;
typedef Kernel::Vector_3                                        Vector;
typedef Kernel::Segment_3                                       Segment;

typedef CGAL::Aff_transformation_3<Kernel>  Aff_transformation;
typedef CGAL::Bbox_3                        Bbox_3;

// Polyhedron
typedef CGAL::Polyhedron_3<Kernel>        Polyhedron;
typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
typedef typename Polyhedron::Facet_handle 	                    Facet_handle;
typedef typename Polyhedron::Facet_iterator                     Facet_iterator;
typedef typename Polyhedron::Vertex_iterator                    Vertex_iterator;
typedef typename Polyhedron::Halfedge_around_facet_circulator   Halfedge_around_facet_circulator;
typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

// basic structures
typedef std::vector<Vertex_handle>          VertexList;
typedef std::set<Vertex_handle>             VertexSet;
typedef std::vector<int>                    IntList;
typedef std::vector<Point>                  PointList;

// KDTree
typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type                         Vertex_point_pmap;
typedef CGAL::Search_traits_3<Kernel>                                                       KDTraits_base;
typedef CGAL::Search_traits_adapter<Vertex_handle, Vertex_point_pmap, KDTraits_base>        KDTraits;
typedef CGAL::Orthogonal_k_neighbor_search<KDTraits>                                        K_neighbor_search;
typedef K_neighbor_search::Tree                                                             KDTree;
typedef KDTree::Splitter                                                                    KDSplitter;

#endif