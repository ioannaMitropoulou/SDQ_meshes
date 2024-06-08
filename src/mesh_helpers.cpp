//
// Created by Ioanna Mitropoulou on 29.11.21.
//

#include "mesh_helpers.h"
#include "helpers.h"

#include <igl/cut_mesh.h>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>
#include <igl/massmatrix.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <directional/cut_mesh_with_singularities.h>
//#include <progresscpp/ProgressBar.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <directional/dcel.h>





SparseMatrix<double> MeshHelpers::VtoF(const MatrixXd& V, const MatrixXi& F) {
    SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    vector<Triplet<double>> data;

    for (int fi = 0; fi < F.rows(); ++fi) {
        double total_mass = 0.0;
        for (int m = 0; m < 3; ++m) {
            int vi = F(fi, m);
            total_mass += M.coeff(vi, vi);
        }
        for (int m = 0; m < 3; ++m) {
            int vi = F(fi, m);
            data.emplace_back(fi, vi, M.coeff(vi, vi) / total_mass);
        }
    }

    SparseMatrix<double> VtoF(F.rows(), V.rows());
    VtoF.setFromTriplets(data.begin(), data.end());
    return VtoF;
}


SparseMatrix<double> MeshHelpers::FtoV(const MatrixXd& V, const MatrixXi& F) {
    VectorXd A;
    igl::doublearea(V, F, A);
    A *= 0.5;

    vector<vector<int>> VF(V.rows());
    for (int vi = 0; vi < V.rows(); ++vi)
        for (int fi = 0; fi < F.rows(); ++fi)
            for (int m = 0; m < 3; ++m)
                if (F(fi, m) == vi)
                    VF[vi].push_back(fi);

    vector<Triplet<double>> data;

    for (int vi = 0; vi < V.rows(); ++vi) {
        double total_vA = 0;
        for (int fi : VF[vi])
            total_vA += A[fi];

        for (int fi : VF[vi])
            data.emplace_back(vi, fi, A[fi] / total_vA);
    }

    SparseMatrix<double> FtoV(V.rows(), F.rows());
    FtoV.setFromTriplets(data.begin(), data.end());
    return FtoV;
}


SparseMatrix<double> MeshHelpers::EtoF(const MatrixXi& F, const MatrixXi& EF) {
    vector<Triplet<double>> data;

    for (int ei=0; ei<EF.rows(); ++ei){
        int f1 = EF(ei,0);
        int f2 = EF(ei, 1);
        if (f1 >=0 && f2 >= 0) {
            data.emplace_back(f1, ei, 1.0);
            data.emplace_back(f2, ei, -1.0);
        }
    }

    SparseMatrix<double> EtoF(F.rows(), EF.rows());
    EtoF.setFromTriplets(data.begin(), data.end());
    return EtoF;
}


int MeshHelpers::edge_index(int v1, int v2, const vector<vector<int>>& VE)
{ // returns the index
    assert(v1 < VE.size());
    assert(v2 < VE.size());
    for (int e1i : VE[v1])
        for (int e2i : VE[v2])
            if (e1i == e2i)
                return e1i;
    return -1;
}


bool MeshHelpers::halfedge_index(int v0, int v1, int ei, const MatrixXi& F, const MatrixXi& EF, int& fi, int& m)
{
    if (ei >= 0) {
        for (int k = 0; k < 2; ++k) {
            fi = EF(ei, k);
            if (fi >= 0) {
                for (m = 0; m < 3; ++m) {
                    if (F(fi, m) == v0 && F(fi, (m + 1) % 3) == v1 ||
                        F(fi, m) == v1 && F(fi, (m + 1) % 3) == v0) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


bool MeshHelpers::halfedge_index(int v0, int v1, const vector<vector<int>>& VE, const MatrixXi& F, const MatrixXi& EF, int& fi, int& m)
{
    const int ei = MeshHelpers::edge_index(v0, v1, VE);

    if (ei >= 0) {
        for (int k = 0; k < 2; ++k) {
            fi = EF(ei, k);
            if (fi >= 0) {
                for (m = 0; m < 3; ++m) {
                    if (F(fi, m) == v0 && F(fi, (m + 1) % 3) == v1 ||
                        F(fi, m) == v1 && F(fi, (m + 1) % 3) == v0) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


void MeshHelpers::complex_local_to_cartesian_world_vf(const MatrixXd& B1, const MatrixXd& B2, MatrixXd& cartesian_vf, const VectorXcd& complex_vf)
{
    int nF = B1.rows();
    cartesian_vf.setZero(nF, 3);
    for (int fi = 0; fi < nF; ++fi)
        cartesian_vf.row(fi) = B1.row(fi) * complex_vf[fi].real() + B2.row(fi) * complex_vf[fi].imag();
}


void MeshHelpers::complex_local_to_cartesian_world_vf(const MatrixXd& B1, const MatrixXd& B2, Vector3d& cartesian_vec, const complex<double>& complex_vec, int fi) {
    cartesian_vec = B1.row(fi) * complex_vec.real() + B2.row(fi) * complex_vec.imag();
}


void MeshHelpers::cartesian_world_to_complex_local_vf(const MatrixXd& B1, const MatrixXd& B2, const MatrixXd& cartesian_vf,  VectorXcd& complex_vf)
{
    int nF = B1.rows();
    complex_vf.setZero(nF);
    for (int fi=0; fi<nF; ++fi)
        complex_vf[fi] = complex<double>(B1.row(fi).dot(cartesian_vf.row(fi)), B2.row(fi).dot(cartesian_vf.row(fi)));
}


void MeshHelpers::cartesian_world_to_complex_local_vf(const MatrixXd& B1, const MatrixXd& B2, const Vector3d& cartesian_vec, complex<double>& complex_vec, int fi)
{
	complex_vec = complex<double>(B1.row(fi).dot(cartesian_vec), B2.row(fi).dot(cartesian_vec));
}


void MeshHelpers::insert_path_to_halfedges_cut_matrix(const vector<int>& path, const MatrixXi& F, const MatrixXi& EF, const vector<vector<int>>& VE, const MatrixXi& TT, 
    const MatrixXi& TTi, MatrixXi& C)
{
    for (int ii = 0; ii < path.size() - 1; ++ii)
    {
        const int v0 = path[ii];
        const int v1 = path[ii + 1];

        int fi, m;
        bool found = MeshHelpers::halfedge_index(v0, v1, VE, F, EF, fi, m);
        assert(found);
        C(fi, m) = 1;
        C(TT(fi, m), TTi(fi, m)) = 1;
    }
}


Vector3d MeshHelpers::barycentric_coordinates(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& p) {
    // Compute barycentric coordinates (u, v, w) for
    // point p with respect to triangle (a, b, c)
    Vector3d v0 = b - a;
    Vector3d v1 = c - a;
    Vector3d v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;

    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;

    return Vector3d(u, v, w);
}


void MeshHelpers::get_cuts_from_singularities(const MatrixXd& V, const MatrixXi& F, const MatrixXi& EF,
    const vector<int>& singularities_vf, const vector<vector<int>>& VE, set<int>& vertices_in_cut, MatrixXi& C) {

    cout << "Cutting mesh from singularities : " << singularities_vf.size() << " singular vertices" << endl;

    std::vector<std::vector<int> > VF, VFi;
    igl::vertex_triangle_adjacency(V, F, VF, VFi);

    std::vector<std::vector<int> > VV;
    igl::adjacency_list(F, VV);

    MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);


    // check how many independent boundaries there are
    vector<vector<Index>> L; // vector of indices of boundary loops
    igl::boundary_loop(F, L); // L: list of loops where L[i] = ordered list of boundary g_vertices in loop i


    for (int i = 0; i < L.size(); ++i)
    {
        if (i == 0) // first boundary: add to the cut all its vertices 
        {
            for (int vi : L[i])
                vertices_in_cut.insert(vi);
        }
        else // other boundaries
        {
            // connect the first vertex of the boundary with the existing cut finding the sortest path from all available shortest paths
            std::vector<int> best_path;
            for (int k = 0; k < L[i].size(); ++k) // try all vertices of the boundary, and select the one that produces the shortest path
            {
                std::vector<int> path;
                VectorXd min_distance;
                VectorXi previous;
                int vertex_found = igl::dijkstra(int(L[i][k]), vertices_in_cut, VV, min_distance, previous);
                igl::dijkstra(vertex_found, previous, path);

                if (path.size() < best_path.size() || best_path.empty())
                    best_path = path;
            }

            vertices_in_cut.insert(best_path.begin(), best_path.end());

            // insert path to cut
            if (!best_path.empty()) insert_path_to_halfedges_cut_matrix(best_path, F, EF, VE, TT, TTi, C);

            // then add all other vertices of boundary on the cut
            for (int vi : L[i]) vertices_in_cut.insert(vi);
        }
    }


    //then, add all singularities one by one by using Dijkstra's algorithm
    for (int i = 0; i < singularities_vf.size(); ++i)
    {
        if (vertices_in_cut.find(singularities_vf[i]) == vertices_in_cut.end()) { // if the singularity is not already in the vertices_in_cut set

			// --- Dijkstra
            std::vector<int> path;
            VectorXd min_distance;
            VectorXi previous;

            if (vertices_in_cut.empty())
                vertices_in_cut.insert(singularities_vf[i]);
            else
            {
                int vertex_found = igl::dijkstra(singularities_vf[i], vertices_in_cut, VV, min_distance, previous);
                igl::dijkstra(vertex_found, previous, path);
            }
            vertices_in_cut.insert(path.begin(), path.end());

            // insert path to cut
            if (!path.empty()) insert_path_to_halfedges_cut_matrix(path, F, EF, VE, TT, TTi, C);
        }
    }
}


void MeshHelpers::get_cuts_edge_correspondence(const MatrixXd& V, const MatrixXi& F, const MatrixXi& EV, const MatrixXi& EF, const VectorXi& vertex_correspondence,
    const vector<vector<int>>& VE, VectorXi& ei_correspondence)
{

    // vertex correspondence 
    vector<vector<int>> cutV_correspondence(V.rows());
    for (int vi = 0; vi < V.rows(); ++vi) {
        for (int vii = vi + 1; vii < V.rows(); ++vii) {
            if (vertex_correspondence[vii] == vertex_correspondence[vi]) {
                cutV_correspondence[vi].push_back(vii);
                cutV_correspondence[vii].push_back(vi);
                break;
            }
        }
    }

    // edge correspondence
    ei_correspondence.setConstant(EV.rows(), -1);
    for (int ei = 0; ei < EV.rows(); ++ei) {
        int v0 = EV(ei, 0);
        int v1 = EV(ei, 1); // current edge vertices

        // if both vertices have a correspondence
        if (!cutV_correspondence[v0].empty() && !cutV_correspondence[v1].empty()) {

            for (int k = 0; k < cutV_correspondence[v0].size(); ++k) { // find other matching edge (ei_other)
                for (int l = 0; l < cutV_correspondence[v1].size(); ++l) {
                    int v0_other = cutV_correspondence[v0][k];
                    int v1_other = cutV_correspondence[v1][l];
                    // find edge between v0_other and v1_other if it exists
                    int ei_other = edge_index(v0_other, v1_other, VE);
                    if (ei_other >= 0) {
                        assert(EF(ei, 0) < 0 || EF(ei, 1) < 0); // make sure ei is a boundary edge
                        assert(EF(ei_other, 0) < 0 || EF(ei_other, 1) < 0); // make sure ei_other is a boundary edge
                        ei_correspondence[ei] = ei_other;
                        ei_correspondence[ei_other] = ei;
                        goto ctn;
                    }
                }
            }
        }

        // if only one of the two vertices have a correspondence
        else if (!cutV_correspondence[v0].empty() || !cutV_correspondence[v1].empty()) {
            int v = cutV_correspondence[v0].empty() ? v0 : v1;
            vector<int> other_vs = cutV_correspondence[v0].empty() ? cutV_correspondence[v1] : cutV_correspondence[v0];
            for (int k = 0; k < other_vs.size(); ++k) {
                int v_other = other_vs[k];
                int ei_other = edge_index(v, v_other, VE);
                if (ei_other >= 0) {
                    assert(EF(ei, 0) < 0 || EF(ei, 1) < 0); // make sure ei is a boundary edge
                    assert(EF(ei_other, 0) < 0 || EF(ei_other, 1) < 0); // make sure ei_other is a boundary edge
                    ei_correspondence[ei] = ei_other;
                    ei_correspondence[ei_other] = ei;
                    goto ctn;
                }

            }
        }
    ctn:;

    }
}


void MeshHelpers::collapse_boundary(MatrixXd& Vc, MatrixXi& Fc, VectorXi vis_to_collapse_group, VectorXi cutV_correspondence,
    VectorXi& F_to_Fc, VectorXi& V_to_Vc)
    // it is intentional that vis_to_collapse_group and cutV_correspondence are not passed by reference; I need to be able to modify them without changing the original
{
    cout << "Collapsing boundary vertices" << endl;
    assert(vis_to_collapse_group.size() == Vc.rows());

    // maps
    const int initial_nV = Vc.rows();
    const int initial_nF = Fc.rows();
    F_to_Fc.setZero(initial_nF);
    V_to_Vc.setZero(initial_nV);
    for (int fi = 0; fi < initial_nF; ++fi) F_to_Fc[fi] = fi;
    for (int vi = 0; vi < initial_nV; ++vi) V_to_Vc[vi] = vi;


    int nV_to_collapse = 0;
    int boundary_max_id = 0;

    for (int vi = 0; vi < Vc.rows(); ++vi)
        if (vis_to_collapse_group[vi] >= 0) {
            ++nV_to_collapse;
            boundary_max_id = max(boundary_max_id, vis_to_collapse_group[vi]);
        }
    const int n_boundaries = boundary_max_id + 1; // number = last index + 1
    cout << "Number of boundaries = " << n_boundaries << endl;

    if (n_boundaries > 0) // if we have boundaries to collapse
    {

        // --- get centers where boundaries should collapse: This center's position doesn't really matter for the computation; it's just nice to have it for visualizing the result as geometry
        vector<int> boundary_vertex_count(n_boundaries, 0);
        vector<RowVector3d> boundary_centers(n_boundaries, RowVector3d(0, 0, 0));
        for (int vi = 0; vi < Vc.rows(); ++vi) {
            int b_id = vis_to_collapse_group[vi];
            if (b_id >= 0) {
                boundary_vertex_count[b_id] += 1;
                boundary_centers[b_id] += Vc.row(vi);
            }
        }
        for (int b_id = 0; b_id < boundary_centers.size(); ++b_id)
            boundary_centers[b_id] /= boundary_vertex_count[b_id];


        // --- collapsing process
        //progresscpp::ProgressBar progressBar(nV_to_collapse, 50, '#', '-');
        bool went_through_entire_mesh = false;

        while (!went_through_entire_mesh)
        {
            // update edge topology
            MatrixXi EVc, FEc, EFc;
            igl::edge_topology(Vc, Fc, EVc, FEc, EFc);

            for (int fi = 0; fi < Fc.rows(); ++fi) {

                for (int m = 0; m < 3; ++m) {
                    int ei = FEc(fi, m);
                    int other_fi = EFc(ei, 0) == fi ? EFc(ei, 1) : EFc(ei, 0); // adjacent face of fi

                    int vi1 = Fc(fi, m); // edge vertices of fi
                    int vi2 = Fc(fi, (m + 1) % 3);

                    if (other_fi < 0) { // if fi is on boundary

                        int bi1 = vis_to_collapse_group[vi1]; // collapsable island index of vertex
                        int bi2 = vis_to_collapse_group[vi2];

                        if (bi1 >= 0 && bi2 >= 0 && bi1 == bi2) // if both vertices belong to the same collapsable boundary
                        {
                            // then collapse face fi to an edge
                            vector<int> vis_to_collapse;
                            int b_id = bi1; // = bi2,  boundary index of current vertices being collapsed

                            vis_to_collapse.push_back(vi1);
                            vis_to_collapse.push_back(vi2);
                            // RowVector3d new_vertex = (Vc.row(vi1) + Vc.row(vi2)) * 0.5;
                            RowVector3d new_vertex = boundary_centers[b_id];


                            // --- update the position of all vertices correspond to vi1 or vi2 before cutting from singularities
                            // (but don't merge them > they need to remain separate vertices)
                            // again: this is only for visualization purposes! The position of the vertices doesn't really matter, only the changes in topology matter
                            for (int v = 0; v < Vc.rows(); ++v)
                                if (cutV_correspondence[v] == cutV_correspondence[vi1] || cutV_correspondence[v] == cutV_correspondence[vi2])
                                    Vc.row(v) = new_vertex;

                            // --- remove old vertices and old per-vertex data
                            sort(vis_to_collapse.begin(), vis_to_collapse.end(), [](int a, int b)->bool { return a > b; });
                            Helpers::removeRow(Vc, vis_to_collapse[0]);
                            Helpers::removeRow(Vc, vis_to_collapse[1]);

                            Helpers::removeRow(vis_to_collapse_group, vis_to_collapse[0]);
                            Helpers::removeRow(vis_to_collapse_group, vis_to_collapse[1]);

                            Helpers::removeRow(cutV_correspondence, vis_to_collapse[0]);
                            Helpers::removeRow(cutV_correspondence, vis_to_collapse[1]);


                            // --- add new vertex and update per-vertex data
                            int old_nV = Vc.rows();
                            int new_vi = old_nV;

                            Vc.conservativeResize(old_nV + 1, 3);
                            Vc.row(new_vi) = new_vertex;

                            vis_to_collapse_group.conservativeResize(old_nV + 1);
                            vis_to_collapse_group[new_vi] = b_id;

                            cutV_correspondence.conservativeResize(old_nV + 1);
                            cutV_correspondence[new_vi] = -b_id - 1; // correspondence of new vertex does not exist, so we have the convention to assign  the -boundary_index -1
                            // hmm, should I here instead assign the old cutV_correspondence? I'm not sure if this value really matters for the further collapsing

                            // --- update V_to_Vc map
                            for (int vi_map = 0; vi_map < initial_nV; ++vi_map) {
                                int v_current = V_to_Vc[vi_map];
                                if (v_current == vi1 || v_current == vi2) { V_to_Vc[vi_map] = new_vi; }
                                else {
                                    if (v_current > vi2) { V_to_Vc[vi_map] -= 1; }
                                    if (v_current > vi1) { V_to_Vc[vi_map] -= 1; }
                                }
                            }

                            // --- remove collapsed face fi
                            Helpers::removeRow(Fc, fi);

                            // --- update all other faces
                            for (int fi_other = 0; fi_other < Fc.rows(); ++fi_other) {
                                for (int mm = 0; mm < 3; ++mm) {
                                    int v_other = Fc(fi_other, mm); // initial unchanged value of vertex

                                    if (v_other == vi1 || v_other == vi2) {
                                        Fc(fi_other, mm) = new_vi;
                                    }
                                    else {
                                        if (v_other > vi2)                    Fc(fi_other, mm) -= 1;
                                        if (v_other > vi1)                    Fc(fi_other, mm) -= 1;
                                    }
                                }
                            }

                            // --- update F_to_Fc map
                            for (int fi_map = 0; fi_map < initial_nF; ++fi_map) {
                                int f_current = F_to_Fc[fi_map];
                                if (f_current == fi) { F_to_Fc[fi_map] = -1; }
                                else if (f_current > fi) { F_to_Fc[fi_map] -= 1; }
                            }

                            //progressBar.display();
                            //++progressBar;

                            goto restart_search; // restart from scratch after the collapse of each edge
                        }
                    }
                }

                if (fi == Fc.rows() - 1)
                    went_through_entire_mesh = true;

            }
        restart_search:;
        }
        //progressBar.done();
    }
}


void MeshHelpers::get_ground_boundary_islands_indices(const MatrixXd& V, const MatrixXi& F, const VectorXi& vs_to_original_boundary,
    const VectorXi& cutV_correspondence, VectorXi& ground_islands)
{
    cout << "Creating ground islands : ";

    // --- get ground vertices
    vector<int> ground_vertices;
    for (int vi = 0; vi < V.rows(); ++vi)
        if (vs_to_original_boundary[vi] >= 0)  // if the vertices belong to a boundary
            if (V(vi, 2) < 1e-3) // if they are on the ground (z=0)
                ground_vertices.push_back(vi);
    cout << "found " << ground_vertices.size() << " ground vertices, ";


    // --- create adjacency graph between them
    vector<vector<int>> VV;
    igl::adjacency_list(F, VV);

    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, int> MyUGraph;
    typedef boost::graph_traits<MyUGraph>::vertex_descriptor Vertex;
    typedef boost::property_map<MyUGraph, boost::vertex_index_t>::type IndexMap;
    typedef boost::graph_traits<MyUGraph>::vertex_iterator vertex_iter;
    MyUGraph G(ground_vertices.size());

    for (int i = 0; i < ground_vertices.size(); ++i) {
        int vi = ground_vertices[i];

        // add edges from mesh adjacency
        for (int adj_vi : VV[vi]) {
            if (Helpers::in_vector(ground_vertices, adj_vi)) {
                int adj_i = Helpers::element_index(ground_vertices, adj_vi);
                boost::add_edge(i, adj_i, G); // add edge between vi and adjacent_vi
            }
        }

        // add edges from cutV_correspondence
        if (i < ground_vertices.size() - 1)
            for (int i_other = i + 1; i_other < ground_vertices.size(); ++i_other)
                if (cutV_correspondence[ground_vertices[i]] == cutV_correspondence[ground_vertices[i_other]])
                    boost::add_edge(i, i_other, G); // add edge between vi and adjacent_vi
    }

    // --- get connected components
    std::vector<int> component(boost::num_vertices(G));
    size_t num_components = boost::connected_components(G, &component[0]);
    cout << num_components << " connected islands" << endl;

    // --- assign values to ground_islands
    ground_islands.setConstant(V.rows(), -1);
    for (int i = 0; i < ground_vertices.size(); ++i) {
        int vi = ground_vertices[i];
        ground_islands[vi] = component[i];
    }
}


Vector3d MeshHelpers::face_gradient(const MatrixXd& V, const MatrixXi& F, int face_index, const Vector3d& face_normal,
    double face_area, const VectorXd& spinning_form) {
    Vector3d vi = V.row(F(face_index, 0));
    Vector3d vj = V.row(F(face_index, 1));
    Vector3d vk = V.row(F(face_index, 2));
    double fi = 0;
    double fj = fi + spinning_form(3 * face_index);
    double fk = fj + spinning_form(3 * face_index + 1);
    RowVector3d g = (fj - fi) * face_normal.cross(vi - vk) / (2 * face_area) +
        (fk - fi) * face_normal.cross(vj - vi) / (2 * face_area);
    return Vector3d(g); // g.normalized()
}


void MeshHelpers::flip_Emap(VectorXi& Emap)
{
    for (int ei = 0; ei < Emap.size(); ++ei)
        if (Emap[ei] >= 0)
            Emap[ei] = static_cast<int>(!static_cast<bool>(Emap[ei]));
}

void MeshHelpers::flip_UVcoords(MatrixXd& UV)
{
    const VectorXd col0 = UV.col(0);
    UV.col(0) = UV.col(1);
    UV.col(1) = col0;
}


void MeshHelpers::get_boundary_polyline_and_vertices(MatrixXd& boundaryPolyline, VectorXi& boundary_vis, const MatrixXd& V, const MatrixXi& F)
{
    int nV = V.rows();
    std::vector<std::vector<Index> > Loops;
    igl::boundary_loop(F, Loops);
    boundary_vis.setZero(nV);
    for (const auto& loop : Loops) {
        for (int vi : loop) {
            boundary_vis[vi] = 1;
        }

        // add boundary polyline segments
        const int existing_boundary_segments = boundaryPolyline.rows();
        boundaryPolyline.conservativeResize(existing_boundary_segments + loop.size(), 6);
        for (int i = 0; i < loop.size() - 1; ++i) {
            boundaryPolyline.block(existing_boundary_segments + i, 0, 1, 3) = V.row(loop[i]);
            boundaryPolyline.block(existing_boundary_segments + i, 3, 1, 3) = V.row(loop[i + 1]);
        }
        boundaryPolyline.block(existing_boundary_segments + loop.size() - 1, 0, 1, 3) = V.row(loop[loop.size() - 1]);
        boundaryPolyline.block(existing_boundary_segments + loop.size() - 1, 3, 1, 3) = V.row(loop[0]);
    }
}


Vector3d MeshHelpers::projectToPolyline(const Vector3d& pt, const MatrixXd& Polyline)
{
    int nB = Polyline.rows();
    VectorXd ds(nB);
    for (int i = 0; i < nB; ++i)
    {
        Vector3d A(Polyline(i, 0), Polyline(i, 1), Polyline(i, 2));
        Vector3d B(Polyline(i, 3), Polyline(i, 4), Polyline(i, 5));
        ds[i] = (pt - Geometry::project_pt_on_line_segment(A, B, pt)).squaredNorm();
    }
    VectorXd::Index min_i;
    ds.minCoeff(&min_i);

    Vector3d A(Polyline(min_i, 0), Polyline(min_i, 1), Polyline(min_i, 2));
    Vector3d B(Polyline(min_i, 3), Polyline(min_i, 4), Polyline(min_i, 5));
    return Geometry::project_pt_on_line_segment(A, B, pt);
}

pair<Vector3d, bool> MeshHelpers::projectToBoundaryInDirection(const Vector3d& pt, const MatrixXd& boundaryPolyline, const Vector3d& dir)
{
    int nB = boundaryPolyline.rows();
    Vector3d D = pt + dir * 100;

    for (int i = 0; i < nB; ++i) {
        Vector3d A(boundaryPolyline(i, 0), boundaryPolyline(i, 1), boundaryPolyline(i, 2));
        Vector3d B(boundaryPolyline(i, 3), boundaryPolyline(i, 4), boundaryPolyline(i, 5));
        double tol = 1e-3;
        pair<Vector3d, bool> pair = Geometry::intersection_line_line(A, B, pt, D, tol);
        if (pair.second)
        {
            return pair;
        }
    }
    return std::make_pair(Vector3d(), false);
}
