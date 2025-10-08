#ifndef QEM_FUNCTION_H
#define QEM_FUNCTION_H

#include "types.h"
#include "metric.h"
#include "pqueue.h"
#include "candidate.h"

#include <map>
#include <string>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/circulator.h>

typedef typename std::map<Facet_handle, QEM_metric>     FacetQemMap;
typedef typename std::map<Vertex_handle, QEM_metric>    VertQemMap;
typedef typename std::map<Vertex_handle, int>           VertIntMap;
typedef typename std::vector<QEM_metric>                QemList;

typedef typename Candidate<Vertex_handle>                   CCandidate;
typedef typename Candidate_more<CCandidate>                 More;
typedef typename Custom_priority_queue<CCandidate, More>    PQueue;

typedef typename QEM_metric::VectorXd                     VectorXd;
typedef typename QEM_metric::MatrixXd                     MatrixXd;


void compute_init_poles(    const Polyhedron& simp_mesh, 
                            const KDTree& tree,
                            VertexSet& init_poles)
{
    init_poles.clear();

    for(Vertex_iterator vd: vertices(simp_mesh))
    {
        Point query = vd->point();
        K_neighbor_search search(tree, query, 1);
        K_neighbor_search::iterator it = search.begin();
        init_poles.insert(it->first);
        // for debug
        //std::cout << "source: " << query << std::endl;
        //std::cout << "target: " << (it->first)->point() << std::endl;
    }
}

void init_pole_colors(int num_poles, IntList& pole_colors)
{
    for(int i = 0; i < num_poles * 3; i++)
    {
        double value = (double) rand() / (RAND_MAX);
        int color = std::max(0, std::min(255, (int)std::round(value * 255)));
        pole_colors.push_back(color);
    }
}

double compute_radius(Polyhedron& mesh)
{
    Bbox_3 box = CGAL::Polygon_mesh_processing::bbox(mesh);
    double x = (box.xmin() + box.xmax()) * 0.5;
    double y = (box.ymin() + box.ymax()) * 0.5;
    double z = (box.zmin() + box.zmax()) * 0.5;
    double r = std::sqrt(CGAL::squared_distance(Point(x, y, z), Point(box.xmin(), box.ymin(), box.zmin())));
    return r;
}

void assemble_qem_map_for_facets(const Polyhedron& orig_mesh, FacetQemMap& facet_qem_map)
{
    for(Facet_handle fd: faces(orig_mesh))
    {
        Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, orig_mesh);
        double area = CGAL::Polygon_mesh_processing::face_area(fd, orig_mesh);
        Halfedge_around_facet_circulator vit = fd->facet_begin();
        Point p = vit->vertex()->point();
        QEM_metric facet_qem;
        facet_qem.init_qem_metrics_face(area, p, normal);
        facet_qem_map.insert({fd, facet_qem});
    }
    std::cout << "Initialized facet qem map size: " << facet_qem_map.size() << std::endl;
}

void assemble_qem_map_for_vertices(const Polyhedron& orig_mesh, FacetQemMap& facet_qem_map, VertQemMap& vert_qem_map)
{
    for(Vertex_iterator vd: vertices(orig_mesh))
    {
        QEM_metric vert_qem;
        Halfedge_around_vertex_circulator vcirc = vd->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;

        CGAL_For_all(vcirc, vend)
        {
            Facet_handle fd = vcirc->facet();
            vert_qem = vert_qem + facet_qem_map[fd];
        }

        vert_qem_map.insert({vd, vert_qem});
    }
    std::cout << "Initialized vertex qem map size: " << vert_qem_map.size() << std::endl;
}


double compute_minimum_qem_error(Point center_point, QEM_metric& query_qem)
{
    VectorXd center_vec(4);
    center_vec << center_point.x(), center_point.y(), center_point.z(), 1.;
    double cost = center_vec.transpose() * query_qem.get_4x4_matrix() * center_vec;
    return std::abs(cost);
}

double compute_collapse_loss(Vertex_handle vd, QEM_metric& cluster_qem, Point& cluster_center, int cluster_label, double dist_weight)
{
    double qem_cost = compute_minimum_qem_error(vd->point(), cluster_qem);
    double dist_cost = CGAL::squared_distance(vd->point(), cluster_center);
    double cost = qem_cost + dist_weight * dist_cost;
    return cost;
}

void region_growing(const VertexSet& init_poles, 
                    VertQemMap& vert_qem_map, 
                    double dist_weight,
                    VertIntMap& vert_labels,
                    QemList& pole_qems,
                    IntList& pole_counts,
                    PointList& pole_centers)
{
    PQueue growing_queue;

    // init poles
    int label = 0;
    for(auto vd: init_poles)
    {
        vert_labels.insert({vd, label});
        pole_qems.push_back(vert_qem_map[vd]);
        pole_counts.push_back(1);
        pole_centers.push_back(vd->point());

        Halfedge_around_vertex_circulator vcirc = vd->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend)
        {
            Vertex_handle oppo_vd = vcirc->opposite()->vertex();
            if(vert_labels.find(oppo_vd) == vert_labels.end())
            {
                double loss = compute_collapse_loss(oppo_vd, pole_qems[label], pole_centers[label], label, dist_weight);
                growing_queue.push(CCandidate(oppo_vd, label, loss));
            }
        }

        label++;
    }

    // growing
    while(!growing_queue.empty())
    {
        CCandidate candidate = growing_queue.top();
        growing_queue.pop();
        Vertex_handle vd = candidate.handle();
        int label = candidate.index();

        // check availability
        if(vert_labels.find(vd) != vert_labels.end())
            continue;
        // assign label
        vert_labels.insert({vd, label});
        // update qem
        pole_qems[label] = pole_qems[label] + vert_qem_map[vd];
        pole_counts[label]++;
        // check neighbor
        Halfedge_around_vertex_circulator vcirc = vd->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend)
        {
            Vertex_handle oppo_vd = vcirc->opposite()->vertex();
            if(vert_labels.find(oppo_vd) == vert_labels.end())
            {
                double loss = compute_collapse_loss(oppo_vd, pole_qems[label], pole_centers[label], label, dist_weight);
                growing_queue.push(CCandidate(oppo_vd, label, loss));
            }
        }
    }

    std::cout << "vert label map size: " << vert_labels.size() << std::endl;
}

Point compute_optimal_point(QEM_metric& cluster_qem, Point& cluster_pole)
{
    // solve Qx = b
    MatrixXd qem_mat = cluster_qem.get_4x4_svd_matrix();
    VectorXd qem_vec = qem_mat.row(3); // 0., 0., 0., 1.
    VectorXd optim(4);

    // check rank
    Eigen::FullPivLU<MatrixXd> lu_decomp(qem_mat);
    lu_decomp.setThreshold(1e-5);

    // full rank -> direct inverse
    if(lu_decomp.isInvertible())
    {
        optim = lu_decomp.inverse() * qem_vec;
    }
    else
    {   // low rank -> svd pseudo-inverse
        Eigen::JacobiSVD<MatrixXd> svd_decomp(qem_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd_decomp.setThreshold(1e-5);

        optim(0) = cluster_pole.x();
        optim(1) = cluster_pole.y();
        optim(2) = cluster_pole.z();
        optim(3) = 1.;

        optim = optim + svd_decomp.solve(qem_vec - qem_mat * optim);
    }

    Point optim_point(optim(0), optim(1), optim(2));

    return optim_point;
}

void compute_optimal_points(PointList& pole_centers, 
                            QemList& pole_qems, 
                            const IntList& pole_counts,
                            PointList& optimal_points)
{
    for(int i = 0; i < pole_qems.size(); i++)
    {
        Point center;

        if(pole_counts[i] == 1)
            center = pole_centers[i];
        else
            center = compute_optimal_point(pole_qems[i], pole_centers[i]);

        optimal_points.push_back(center);
    }
}

void update_poles(const PointList& new_pole_centers,
                  const KDTree& tree,
                  VertexSet& new_poles)
{
    for(int i = 0; i < new_pole_centers.size(); i++)
    {
        Point query = new_pole_centers[i];
        K_neighbor_search search(tree, query, 1);
        K_neighbor_search::iterator it = search.begin();
        new_poles.insert(it->first);
        // for debug
        //std::cout << "source: " << query << std::endl;
        //std::cout << "target: " << (it->first)->point() << std::endl;
    }
}

void save_poles(VertexSet& poles, std::string filename)
{
    PointList pole_points;
    for(auto vd: poles)
        pole_points.push_back(vd->point());
    CGAL::IO::write_XYZ(filename, pole_points);
}

void save_cluster_color_ply(Polyhedron& mesh, VertIntMap& vert_labels, IntList& pole_colors, std::string filename)
{
    std::ofstream cluster_file;
    cluster_file.open(filename);

    cluster_file << "ply\n"
                << "format ascii 1.0\n"
                << "element vertex " << mesh.size_of_vertices() << "\n"
                << "property float x\n"
                << "property float y\n"
                << "property float z\n"
                << "property uchar red\n"
                << "property uchar green\n"
                << "property uchar blue\n"
                << "end_header\n";

    for(Vertex_handle vd: vertices(mesh))
    {
        Point p = vd->point();

        int label = -1;
        if(vert_labels.find(vd) != vert_labels.end())
            label = vert_labels[vd];

        cluster_file << p.x() << " " << p.y() << " " << p.z() << " ";

        if(label == -1)
            cluster_file << "0 0 0" << std::endl;
        else
            cluster_file  << pole_colors[label * 3] << " " << pole_colors[label * 3 + 1] << " " << pole_colors[label * 3 + 2] << std::endl;
    }

    cluster_file.close();
}


#endif