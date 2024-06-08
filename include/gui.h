//
//  Created by Ioanna Mitropoulou on 14.09.22.
//

#pragma once
#ifndef NP_3DP_GUI_H
#define NP_3DP_GUI_H

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>


#include "Np3dpContext.hpp"
#include "vector_field.h"

using namespace std;
using namespace Eigen;

class Gui
{
public:
	Gui(const shared_ptr<Np3dpContext>& np3dp, const shared_ptr<VectorField>& vf): np3dp(np3dp), vf(vf)
	{
		viewer.callback_key_down = [this](auto& viewer, unsigned char key, int modifier) { return key_down(viewer, key, modifier); };
		viewer.callback_mouse_down = [this](auto& viewer, unsigned char key, int modifier) { return mouse_down(viewer, key, modifier); };

		set_mesh_visual(np3dp->mesh.V, np3dp->mesh.F);
		viewer.data().shininess = 10e7;
		update_mesh_colors();
		display_vector_field();

		viewer.plugins.push_back(&plugin);
		plugin.widgets.push_back(&menu);


		// Add menu boxes
		draw_menu_boxes();

		// Launch the viewer
		bool fullscreen = false;
		int windowWidth = 2400;// 1400;
		int windowHeight = 1200;//800;
		viewer.launch(fullscreen, "np3dp", windowWidth, windowHeight);
	}

	bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
	bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);

	void curl_reduction_one_step();


	void set_mesh_visual(const MatrixXd& V, const MatrixXi& F);
	void update_mesh_colors();
	void display_face_degree_mesh_colors();
	void display_fvf_constraints_face_colors();
	void prepare_viewer_index_data(int& index);
	void display_per_vertex_value(const VectorXd& values);
	void display_Emap(const MatrixXd& Vquad, const MatrixXi& EVquad, const VectorXi& Emap);
	void display_edge_maps();
	void mark_edges(const MatrixXd& V, const MatrixXi& EV, const vector<int>& eis, int& viewer_index, const RowVector3d& color, double line_width = 10.0);
	void mark_edges(int& viewer_index, const RowVector3d& color, const MatrixXd& P1, const MatrixXd& P2, double line_width = 10.0);
	void mark_points(int& viewer_index, const RowVector3d& color, const MatrixXd& P);
	// void mark_fixed_vis();
	void mark_quad_edges();
	void display_mesh_boundaries(const MatrixXd& V, const MatrixXi& F);
	void display_vf_singularities();
	void display_quad_mesh_vertex_degrees();
	void display_vector_field();
	void display_paths();
	void display_UV_coords(const MatrixXd& VMeshCut, const MatrixXi& FMeshCut, const MatrixXd& UVcoords);
	void display_partitioned_mesh();
	void draw_partitioning_cuts();
	void display_selected_strip_data();
	void display_strips_centerlines(int& viewer_index, const VectorXi& strip, const RowVector3d& color, int direction);
	void display_markings(const vector<Vector3d>& positions);
	void display_selected_vis(const RowVector3d& color = RowVector3d(1.0, 0.0, 0.0));
	void display_singularity_alignments(int& viewer_index, const MatrixXd& P1, const MatrixXd& P2, const vector<int>& vis, const RowVector3d& color = RowVector3d(0.15, 0.10, 0.15), const RowVector3d& other_color = RowVector3d(0.75, 0.70, 0.75));
	void display_blocked_strips();
	void display_intermediate_strips_of_selected_vis(int& viewer_index);
	void update_visualization_after_topology_edits();


	// -----------------------> Evaluation of quad mesh results
    void evaluate_edges_length_uniformity(int dir = 0);
    void evaluate_orthogonality();
    void evaluate_alignment();
    VectorXd measure_quad_planarity_error(bool set_non_quads_to_zero);

	template<class V>
	inline void display_per_face_data_on_quad_mesh(const MatrixXi& Tquad, const VectorXi& TF, const V& Fmap, bool normalize=true)
	{
		int nF = Tquad.rows();
	    assert(TF.rows() == nF);
	    MatrixXd face_colors(nF, 3);

	    V Map_triangles(nF);
	    for (int fi = 0; fi < nF; ++fi)
	        Map_triangles[fi] = Fmap[TF[fi]];

	    igl::parula(Map_triangles, normalize, face_colors);
	    for (int fi = 0; fi < nF; ++fi)
	        if (Map_triangles[fi] == 0)
	            face_colors.row(fi) = np3dp->mesh_base_color;
	    viewer.data().set_colors(face_colors);
	}

	void serialize() const; 
	void deserialize();
	void invalidate(bool keep_vf_params = false);


	void draw_menu_boxes();

private:
	shared_ptr<Np3dpContext> np3dp;
	shared_ptr<VectorField> vf;
	
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiPlugin plugin;
	igl::opengl::glfw::imgui::ImGuiMenu menu;

	// --- viewer indices
	int vf_viewer_index = -1;
	int vf_rotated_viewer_index = -1;
	int vf_sings_viewer_index = -1;
	int vf_streamlines_viewer_index = -1;
	int mesh_boundaries_viewer_index = -1;
	int paths_viewer_index = -1;
	int strip_sequences_viewer_index = -1;
	int quad_edges_viewer_index = -1;
	int quad_vertex_degrees_viewer_index = -1;
	int emap_viewer_index = -1;
	int cuts_on_edges_viewer_index = -1;
	int mark_selected_viewer_index = -1;
	int mark_selected_vertices_viewer_index = -1;
	int mark_fixed_pts_viewer_index = -1;
	int mark_helper_strips_viewer_index = -1;
	int alignment_of_selected_vertices_viewer_index = -1;
	int position_markers_viewer_index = -1;
	int piece_labels_viewer_index = -1;

	// --- display parameters
	bool show_mesh_boundaries = false;
	bool show_vf_sings = true;
	bool show_vector_field = true;
	bool show_paths = true;
	bool show_paths_rotated = true;
	bool show_quad_vertex_degrees = false;
	bool show_helper_strips = true;
	bool show_helper_counter_strips = true;
	bool show_edge_maps = true;
	bool show_ribs = false;
	bool show_vfConstraints_mesh_colors = false;
	bool show_face_degree_mesh_colors = true;
	bool show_quad_partition_map = false;
	bool show_singularities_alignment = true;

	// --- gui selection parameters
	bool strip_selection_is_on = false;
	bool vertex_selection_is_on = false;
	bool edge_selection_is_on = false;
	bool face_selection_is_on = false;

	bool mark_edge_loop = false; // when selecting an edge, if the entire edge sequence should be marked
	bool save_graphviz = false; // save graphviz visualization of strip graphs. Resulting .png image saved in the output folder.

	// --- mesh colors
	enum MeshColors { None = 0, VfConstraints, FaceDegree, QuadPartitionMap };
	MeshColors meshColors = VfConstraints;
};

#endif //NP_3DP_GUI_H