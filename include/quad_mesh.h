//
// Created by Ioanna Mitropoulou on 07.01.22.
//

#pragma once
#ifndef NP_3DP_QUAD_MESH_H
#define NP_3DP_QUAD_MESH_H

#include <directional/setup_integration.h>
#include <directional/TriMesh.h>
#include <Eigen/Core>
#include "mesh_helpers.h"
#include "structs.h"

using namespace Eigen;
using namespace std;


class QuadMesh
{
public:
    QuadMesh(){}
    QuadMesh(const std::string& output_folder) : output_folder(output_folder){}
    void save_data() const;

    MatrixXd V;
	MatrixXi F, T; // Fquad: quad faces, Tquad: quad triangulated faces 
	VectorXi D, TF; // D: #Fx1 vector that stores the degree (amount of vertices or edges) in a face, TF: relates every triangle to the face it tesselates

	VectorXi Emap; // always for strip_network=0. (Equad x 1): stores 0/1 for each edge if it is aligned with the dominant/subdominant direction respectively. all boundary edges have Emap = -1 
    MatrixXd UVcoords; // always for strip_network=0 (Vquad x 2)

    VectorXi ERibMap; // #Equad x 1: stores 0/1 for each edge that is/is not a rib

    MatrixXd faceNormals;

    // edge topology data
    MatrixXi EV, FE, EF, EFi;
	MatrixXd FEs;
    VectorXi innerEdges;
    VectorXi vertexDegree;

    // boundaries
    VectorXi boundaryVertices; // Vquad x 1
    VectorXi boundaryFaces; // Fquad x 1 : 0/1 if face isn't/is on the boundary

    // halfedge data
    MatrixXi EH, FH;
	VectorXi VH, HE, HF, HV, nextH, prevH, twinH;

    // folders
    string output_folder;

    void generate_from_integration(const directional::TriMesh& mesh, const directional::TriMesh& meshCut, const directional::IntegrationData& intData, const MatrixXd& UVcoords);

    void update_quad_mesh_data();
    bool has_updated_mesh_data() const;

    
    void create_quad_Emap(bool remember_existing_emap_values = false);
    void add_extra_vertex_on_triangle_faces();
    void get_UVcoords_on_quad_vertices(const directional::TriMesh& meshCut, const MatrixXd& triangle_mesh_UVcoords);

	int find_distance_to_boundary(int he) const;
    double get_avg_len_edge_type(int edge_type) const; // returns average length of edges with Emap == edge_type
    MatrixXd get_lines_matrix_from_cuts(const VectorXi& C) const;
    vector<vector<int>> adjacency_list() const;
    vector<vector<int>> adjacency_list_with_Emap(int dir) const;
    bool shortest_path(int source_vertex, const set<int>& target_vertices, const vector<vector<int>>& adjacency_list, vector<int>& path) const;

	static int count_Emap_values_of_face_edges(int fi, const VectorXi& D, const MatrixXi& FE, const VectorXi& Emap, int emap_value_to_count);
	static bool is_boundary_edge(int ei, const MatrixXi& EH);
    static bool is_boundary_face(int f0, const VectorXi& twinH, const VectorXi& D, const MatrixXi& FH);
    static MatrixXd per_face_normals(const MatrixXd& V, const MatrixXi& T, const MatrixXi& F, const VectorXi& TF);
	static VectorXi vertices_on_boundaries(const MatrixXd& V, const MatrixXi& EF, const MatrixXi& EV);
    static VectorXi faces_on_boundaries(const MatrixXi& F, const MatrixXi& EF);
    static VectorXi mesh_vertex_degrees(const MatrixXd& V, const MatrixXi& F, const VectorXi& D);
	static VectorXd face_areas(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const VectorXi& TF);
	static void expand_cutting_from_he(int he, VectorXi& C, vector<int>& passed_vis, int nE, const VectorXi& HE, const VectorXi& HV, const VectorXi& nextH, const VectorXi& twinH, const VectorXi& vertexDegree, const VectorXi& boundaryVertices, bool stop_at_cut_intersections, bool stop_at_irregular_vertices);


	// ------------------ Ribs on edges
	int ribs_dist = 4; // the distance between two consecutive ribs (measured in quads)
	int ribs_f0 = 0; // arbitrary selection of first face where the search for ribs starts from
    void create_quad_ERibs();
    MatrixXd get_ribs_polylines() const;


	// ------------------ manual selection data (vertices, faces, edges, strips)
    SelectionData selection;


    // ------------------ Strip graphs and data for aligning vertices
    StripNetworks strip_networks; 

	vector<pair<int, int>> vis_to_check_for_alignments;
	vector<AlignedVertices> aligned_vertices; // stores aligned singularities

    


    // ------------------ cleanup
    // faces removal: return removed_fis, sorted from highest to lowest. Also remove them from F and D
    vector<int> remove_tiny_faces();
    vector<int> remove_boundary_isolated_triangles(); // removes boundary triangles that have 2 edges on the boundary
    vector<int> remove_faces_without_dominant_edges(); // returns removed_fis

    // vertices removal: return removed_vis, sorted from hightest to lowest. Also remove them from V, UV, FixedVis
    vector<int> cleanup_unreferenced_vertices();
    vector<int> simplify_boundary_polyline();

    void cleanup(); // cleans up the quad mesh


    // ------------------ edit operations
    // Collapse
    void collapse_strip(VectorXi& strip, int selected_strip_direction);
    void collapse_edge(int ei);

    // Smoothing
	bool smoothing_with_boundaries = false;
	int smoothing_iterations = 10;
    void smoothen(const directional::TriMesh& tri_mesh, bool with_smoothened_boundary, int iterations = 10, bool with_cleanup = true); // exceptions_to_fixed_vis: fixed vis which shoudl however be smoothened as free here
    void smoothen_specific_quad_mesh(MatrixXd& V, const MatrixXi& F, const VectorXi& D, const MatrixXd& Triangle_V, const MatrixXi& Triangle_F, bool with_smoothened_boundary, int iterations);

	// Rewiring
	bool rewire_strip_forward = true;
    void rewire_strip(VectorXi& strip, const array<vector<int>, 2>& strip_vis, bool strip_is_closed, bool rewire_forward);
private: // rewiring helpers
    bool find_alignment_of_two_vertex_groups(const VectorXi& strip, const vector<int>& visA, const vector<int>& visB, int& i, int& j, int& fi_common);
    void get_shortest_path_to_vis(int source, vector<int> target_vis, vector<int> other_vis, const vector<vector<int>>& VV, vector<int>& path, bool& path_contains_strip);
    bool decide_rewiring_direction(const vector<int>& visA, const vector<int>& visB);
public:

    // Subdivide
    void subdivide_strip(const VectorXi& strip, int direction, bool save_results=true);
	void subdivide_strip(MatrixXd& V, MatrixXi& F, VectorXi& D, VectorXi& Emap, MatrixXd& UV, const vector<VectorXi*>& F_based_data,
						 const VectorXi& strip, int direction, int& e_remember_v1, int& e_remember_v2, int& e_remember_Emap);
	void subdivide_whole_quad_mesh(bool only_elongated_strips);


    // ------------------ (de)serialize and invalidate
    void serialize(const std::string& serialize_folder) const;
    void deserialize(const std::string& serialize_folder);


	// ------------------ Save prev information to be able to undo a manual change (at the level of not-partitioned quad mesh)
    bool has_previous_data = false;
    void remember_previous_data();
    void restore_previous_data();
    void invalidate_previous_data();
private: // stored previous data
    MatrixXd V_prev;
    MatrixXi F_prev;
    VectorXi D_prev;
    MatrixXd UV_prev;
    VectorXi Emap_prev;
};

#endif //NP_3DP_PARTITIONING_Hboundary
