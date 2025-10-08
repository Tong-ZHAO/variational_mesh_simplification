#include "function.h"

int main(int argc, char* argv[])
{
    if(argc != 4)
	{
		std::cout << "Usage: ./main orig_mesh.obj simplified_mesh.obj dist_weight" << '\n';
		return 1;
	}

    std::string orig_mesh_name = argv[1];
    std::string simp_mesh_name = argv[2];
    double dist_weight = std::stod(argv[3]);

    // read meshes
    Polyhedron orig_mesh, simp_mesh;
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(orig_mesh_name, orig_mesh))
    {
        std::cout << "Error: cannot read file " << orig_mesh_name << "." << std::endl;
        return 1;
    }
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(simp_mesh_name, simp_mesh))
    {
        std::cout << "Error: cannot read file " << simp_mesh_name << "." << std::endl;
        return 1;
    }
    std::cout << "Meshes loaded!" << std::endl;

    std::cout << "Original mesh: " << orig_mesh.size_of_vertices() << " vertices and " << orig_mesh.size_of_facets() << " faces." << std::endl;
    std::cout << "Simplified mesh: " << simp_mesh.size_of_vertices() << " vertices and " << simp_mesh.size_of_facets() << " faces." << std::endl;

    // 0 - initialize kdtree for orig mesh
    Vertex_point_pmap vppmap = get(CGAL::vertex_point, orig_mesh);
    KDTree tree(vertices(orig_mesh).begin(), vertices(orig_mesh).end(), KDSplitter(), KDTraits(vppmap));

    // 0 - change relative dist weight
    double radius = compute_radius(orig_mesh);
    dist_weight = dist_weight * radius;

    // 0.1 - find initial cluster centers for orig mesh
    VertexSet curr_poles;
    compute_init_poles(simp_mesh, tree, curr_poles);
    std::cout << "Initialized " << curr_poles.size() << " clusters!" << std::endl;

    IntList pole_colors;
    init_pole_colors((int)curr_poles.size(), pole_colors);

    // 0.2 - initialize qem map for orig mesh
    FacetQemMap facet_qem_map;
    assemble_qem_map_for_facets(orig_mesh, facet_qem_map);

    VertQemMap vert_qem_map;
    assemble_qem_map_for_vertices(orig_mesh, facet_qem_map, vert_qem_map);

    // clustering
    VertIntMap vert_labels; // map (vertex_handle, clsuter_label)
    QemList pole_qems;      // vector (qem for each cluster)
    IntList pole_counts;    // vector (vertex count for each cluster)
    PointList pole_centers; // vector (center for each cluster)
    PointList optimal_points; // vector (updated center for each cluster)

    for(int i = 0; i < 10; i++)
    {
        std::cout << "Begin " << i+1 << " iteration!" << std::endl;
        vert_labels.clear();
        pole_qems.clear();
        pole_counts.clear();
        pole_centers.clear();
        region_growing(curr_poles, vert_qem_map, dist_weight, 
                       vert_labels, pole_qems, pole_counts, pole_centers);
        
        optimal_points.clear();
        compute_optimal_points(pole_centers, pole_qems, pole_counts, 
                               optimal_points );
        
        curr_poles.clear();
        update_poles(optimal_points, tree,
                     curr_poles);
    }

    // output
    save_poles(curr_poles, "new_poles.xyz");
    save_cluster_color_ply(orig_mesh, vert_labels, pole_colors, "cluster.ply");

    return 0;
}