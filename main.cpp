#include "function.h"

int main(int argc, char* argv[])
{
    if(argc != 6)
	{
		std::cout << "Usage: ./main orig_mesh.obj simplified_mesh.obj dist_weight num_add num_iter" << '\n';
		return 1;
	}

    std::string orig_mesh_name = argv[1];
    std::string simp_mesh_name = argv[2];
    double dist_weight = std::stod(argv[3]);
    int num_add = std::stoi(argv[4]);
    int num_iter = std::stoi(argv[5]);

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

    // 0.2 - initialize qem map for orig mesh
    FacetQemMap facet_qem_map;
    assemble_qem_map_for_facets(orig_mesh, facet_qem_map);

    VertQemMap vert_qem_map;
    assemble_qem_map_for_vertices(orig_mesh, facet_qem_map, vert_qem_map);

    // begin iteration for clustering
    VertIntMap vert_labels; // map (vertex_handle, clsuter_label)
    QemList pole_qems;      // vector (qem for each cluster)
    IntList pole_counts;    // vector (vertex count for each cluster)
    PointList pole_centers; // vector (center for each cluster)
    PointList optimal_points; // vector (updated center for each cluster)

    PairList pole_loss; // vector (loss for each cluster)

    for(int i = 0; i < num_iter; i++)
    {
        std::cout << "Begin " << i+1 << " iteration!" << std::endl;

        for(int j = 0; j < 10; j++)
        {
            std::cout << "  Begin " << j+1 << " clustering!" << std::endl;
            vert_labels.clear();
            pole_qems.clear();
            pole_counts.clear();
            pole_centers.clear();
            region_growing(curr_poles, vert_qem_map, dist_weight, 
                        vert_labels, pole_qems, pole_counts, pole_centers);
            
            if(j != 9)
            {
                optimal_points.clear();
                compute_optimal_points(pole_centers, pole_qems, pole_counts, 
                                    optimal_points );
                
                curr_poles.clear();
                update_poles(optimal_points, tree,
                            curr_poles);
            }
            else if(i != num_iter-1)
            {
                // compute loss for each cluster
                pole_loss.clear();
                compute_cluster_loss(pole_centers, vert_qem_map, vert_labels, 
                                    pole_loss);
                // split cluster with adding num_add poles
                split_cluster(pole_loss, num_add, 
                            curr_poles);
                std::cout << "  Current number of poles: " << curr_poles.size() << std::endl;
            }
        }

    }

    // construct triangle soup
    IntIntList poly_list;
    reconstruct_polygon(orig_mesh, vert_labels, poly_list);
    save_polygon_soup(pole_centers, poly_list, "polygon_soup.off");

    // output
    IntList pole_colors;
    init_pole_colors((int)curr_poles.size(), pole_colors);

    save_poles(curr_poles, "new_poles.xyz");
    save_cluster_color_ply(orig_mesh, vert_labels, pole_colors, "cluster.ply");

    return 0;
}