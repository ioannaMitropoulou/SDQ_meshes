#include "gui.h"
#include "igl/unproject_onto_mesh.h"
#include "quad_mesh.h"
#include <hedra/polygonal_write_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <igl/parula.h>
#include <directional/polygonal_edge_topology.h>
#include <directional/Deprecated/representative_to_raw.h>
#include <igl/point_mesh_squared_distance.h>
#include <boost/filesystem.hpp>
#include <hedra/planarity.h>

#include "imgui.h"
#include "partitioning.h"
#include "paths_tracing.h"
#include "strips.h"


bool Gui::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    if (button == (int)igl::opengl::glfw::Viewer::MouseButton::Right)
        return false;

    if ((strip_selection_is_on || vertex_selection_is_on || edge_selection_is_on || face_selection_is_on) && np3dp->has_quad_mesh)
    {
		SelectionData& sel = np3dp->quad_mesh.selection;

        int triangle_fid;
        Vector3f bc;

        const MatrixXd& V = np3dp->quad_mesh.V;
        const MatrixXi& F = np3dp->quad_mesh.T;

        // Cast a ray in the view direction starting from the mouse position
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, triangle_fid, bc)) {

            if (vertex_selection_is_on)
            {
                int max;
                bc.maxCoeff(&max);
                int selected_vi = F(triangle_fid, max);
                cout << "selected vertex index = " << selected_vi << endl;


                // --- save selected vertex 
                if (sel.selected_vis.size() < sel.n_vis_to_select) {
                    int n = sel.selected_vis.size();
                    sel.selected_vis.resize( n+1);
                    sel.selected_vis[n] = selected_vi;
                }
                else {
                    for (int i=0; i<sel.n_vis_to_select-1; ++i) // push everything one position backwards
                    {
                        sel.selected_vis[i] = sel.selected_vis[i + 1];
                    }
                    sel.selected_vis[sel.n_vis_to_select-1] = selected_vi; // the last slot takes the new value 
                }
                cout << "selected vis : ";  for (int vi : sel.selected_vis) cout << vi << " "; cout << endl;

                // --- display
                if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
                if (alignment_of_selected_vertices_viewer_index >= 0) viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
                prepare_viewer_index_data(mark_selected_vertices_viewer_index);
                viewer.data_list[mark_selected_vertices_viewer_index].point_size = 15; // points not displaying on Windows 11...
                viewer.data_list[mark_selected_vertices_viewer_index].label_size = 5;

				viewer.data_list[mark_selected_vertices_viewer_index].clear_labels();
                for (int k=0; k<sel.selected_vis.size(); ++k)
                {
					int vi = sel.selected_vis[k];
                    viewer.data_list[mark_selected_vertices_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), RowVector3d(1.0, 0.0, 0.0));
	                viewer.data_list[mark_selected_vertices_viewer_index].add_label(np3dp->quad_mesh.V.row(vi), std::to_string(k));
                }
					
                return true;
            }

            if (strip_selection_is_on)
            {
				int& fi = sel.selected_fi;
                fi = np3dp->quad_mesh.TF[triangle_fid];
                cout << "selected face index = " << fi << endl << "face vertices : "; for (int k = 0; k < np3dp->quad_mesh.D[fi]; ++k) cout << np3dp->quad_mesh.F(fi, k) << " "; cout << endl;
            	if (Strips::get_strip_from_fi(np3dp->quad_mesh.D, np3dp->quad_mesh.F, np3dp->quad_mesh.boundaryFaces, np3dp->quad_mesh.Emap, fi, sel.selected_strip_direction, sel.selected_strip, sel.selected_strip_vis, sel.strip_is_closed, vector<int>(), sel.selected_strip_fis)) {
                    cout << "Traced strip with " << sel.selected_strip.sum() << " faces. Strip is closed : " << sel.strip_is_closed << endl;
                    cout << "Press SPACE to flip the strip direction" << endl;
                    if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
                    display_selected_strip_data();
                }
                return true;
            }

            if (face_selection_is_on)
            {
				int& fi = sel.selected_fi;
                fi = np3dp->quad_mesh.TF[triangle_fid];
                cout << "selected face index = " << fi << endl << "face vertices : "; for (int k = 0; k < np3dp->quad_mesh.D[fi]; ++k) cout << np3dp->quad_mesh.F(fi, k) << " "; cout << endl;
                VectorXi Fmap; Fmap.setZero(np3dp->quad_mesh.F.rows()); Fmap[fi] = 1;
                display_per_face_data_on_quad_mesh(np3dp->quad_mesh.T, np3dp->quad_mesh.TF, Fmap); // color mesh faces from strip
                return true;
            }

            if (edge_selection_is_on)
            {
				int& ei = sel.selected_ei;
                ei = -1;

                int max;
                bc.maxCoeff(&max);
                bc[max] = 0.0;
                int second_max;
                bc.maxCoeff(&second_max);
                int v1 = F(triangle_fid, max);
                int v2 = F(triangle_fid, second_max);

                // find edge that connects v1, v2
                for (int k = 0; k < np3dp->quad_mesh.EV.rows(); ++k) // iterate through all edges
                    if (np3dp->quad_mesh.EV(k, 0) == v1 && np3dp->quad_mesh.EV(k, 1) == v2 ||
                        np3dp->quad_mesh.EV(k, 1) == v1 && np3dp->quad_mesh.EV(k, 0) == v2)
                    {
                        ei = k;
                        std::cout << "found edge = " << ei << " with v1 = " << v1 << " , v2 = " << v2 << "and Emap value : " << np3dp->quad_mesh.Emap[ei] << endl;
                        mark_edges(np3dp->quad_mesh.V, np3dp->quad_mesh.EV, { ei }, mark_selected_viewer_index, RowVector3d(1.0, 0.0, 0.0), 12);
                        break;
                    }

                if (ei >= 0 && mark_edge_loop && np3dp->quad_mesh.Emap[ei] >= 0) // mark points of edge loop
                {
                    // get all vertices that belong to the extended edge
                    VectorXi C; C.setZero(np3dp->quad_mesh.EV.rows()); // cuts (#E x 1) 0: not a cut, 1: cut

                    vector<int> passed_vis1, passed_vis2;
					const QuadMesh& q = np3dp->quad_mesh;
                    int he = q.EH(ei, 0);
                    QuadMesh::expand_cutting_from_he(he, C, passed_vis1, C.rows(), q.HE, q.HV, q.nextH, q.twinH, q.vertexDegree, q.boundaryFaces, false, true);
                    he = np3dp->quad_mesh.nextH[np3dp->quad_mesh.twinH[np3dp->quad_mesh.nextH[np3dp->quad_mesh.twinH[he]]]];
                    QuadMesh::expand_cutting_from_he(he, C, passed_vis2, C.rows(), q.HE, q.HV, q.nextH, q.twinH, q.vertexDegree, q.boundaryFaces, false, true);

                    vector<int> vis = passed_vis1;
                    for (int vi : passed_vis2)
                        if (!Helpers::in_vector(vis, vi))
                            vis.push_back(vi);

                    MatrixXd P(vis.size(), 3);
                    int count = 0;
                    for (int vi : vis) {
                        P.row(count) = np3dp->quad_mesh.V.row(vi);
                        ++count;
                    }
                    viewer.data_list[mark_selected_viewer_index].point_size = 10;
                    viewer.data_list[mark_selected_viewer_index].add_points(P, RowVector3d(0.85, 0.15, 0.45));

					viewer.data_list[mark_selected_viewer_index].label_size = 3;
	                for (int k=0; k<vis.size(); ++k)
	                {
		                viewer.data_list[mark_selected_viewer_index].add_label(np3dp->quad_mesh.V.row(vis[k]), std::to_string(k));
	                }
                }
                return true;
            }
        }
        return false;
    }
    return false;
}


bool Gui::key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    cout << "Pressed key : " << key << endl;
    if (key == 's' || key == 'S') {
        //vf->setup();
        //if (show_vector_field) display_vector_field();
        //update_mesh_colors();

        strip_selection_is_on = !strip_selection_is_on;
        if (strip_selection_is_on) {
            vertex_selection_is_on = false; edge_selection_is_on = false; face_selection_is_on = false;
            np3dp->quad_mesh.selection.selected_fi = -1;  np3dp->quad_mesh.selection.selected_ei = -1;
        }
        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
        if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
        if (alignment_of_selected_vertices_viewer_index >= 0) viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
        viewer.data().set_colors(np3dp->mesh_base_color);

        return true;
    }
    if (key == 'v' || key == 'V') {
        if (!vf->is_setup) vf->setup();
        vf->minimize_energy();
        if (show_vector_field) display_vector_field();
        update_mesh_colors();
        return true;
    }
    if (key == 'd' || key == 'D') { // one iteration step
        if (!vf->is_setup) vf->setup();
        vf->direct_solve();
        if (show_vector_field) display_vector_field();
        update_mesh_colors();
        return true;
    }
    if (key == 'c' || key == 'C') { // curl reduction, one step
        if (vf->is_setup)
            curl_reduction_one_step();
        return true;
    }
    if (key == 'q' || key == 'Q') {
        deserialize();
        return true;
    }
    if (key == 'z' || key == 'Z') {
        if (np3dp->quad_mesh.has_previous_data) {
			np3dp->quad_mesh.restore_previous_data();
			np3dp->quad_mesh.strip_networks.invalidate();
            update_visualization_after_topology_edits();
        }
    }
    if (key == 'e' || key == 'E') // collapse strip
    {
        if (strip_selection_is_on) {
		    if (np3dp->quad_mesh.selection.selected_strip.rows() == np3dp->quad_mesh.F.rows()) {
		        np3dp->quad_mesh.collapse_strip(np3dp->quad_mesh.selection.selected_strip, np3dp->quad_mesh.selection.selected_strip_direction);
		        np3dp->quad_mesh.strip_networks.invalidate();
                update_visualization_after_topology_edits();
		    }
        }

        if (edge_selection_is_on){
        	if (np3dp->quad_mesh.selection.selected_ei >= 0) {
                np3dp->quad_mesh.collapse_edge(np3dp->quad_mesh.selection.selected_ei);

                np3dp->quad_mesh.strip_networks.invalidate();

                // update mesh data
                np3dp->quad_mesh.update_quad_mesh_data();
                np3dp->quad_mesh.create_quad_Emap();

                // saving new quad mesh
                igl::writeOBJ(DATA_PATH + np3dp->output_folder + "quad_mesh.obj", np3dp->quad_mesh.V, np3dp->quad_mesh.F); cout << "Saved : " << DATA_PATH + np3dp->output_folder + "quad_mesh.obj" << endl;
                hedra::polygonal_write_OFF(DATA_PATH + np3dp->output_folder + "quad_mesh.off", np3dp->quad_mesh.V, np3dp->quad_mesh.D, np3dp->quad_mesh.F); cout << "Saved : " << DATA_PATH + np3dp->output_folder + "quad_mesh.off" << endl;
                Helpers::write_Emap_to_txt(DATA_PATH + np3dp->output_folder + "Emap.txt", np3dp->quad_mesh.Emap, np3dp->quad_mesh.EV);

                update_visualization_after_topology_edits();
            }
        }
        return true;
    }
    if (key == 'r' || key == 'R') // subdivide strip
    {
        if (np3dp->quad_mesh.selection.selected_strip.rows() == np3dp->quad_mesh.F.rows()) {
            np3dp->quad_mesh.subdivide_strip(np3dp->quad_mesh.selection.selected_strip, np3dp->quad_mesh.selection.selected_strip_direction, true);
	        np3dp->quad_mesh.strip_networks.invalidate();
            update_visualization_after_topology_edits();
        }
        return true;
    }
    if (key == ' ') {
        if (strip_selection_is_on && np3dp->has_quad_mesh) {
            np3dp->quad_mesh.selection.selected_strip_direction = int(!bool(np3dp->quad_mesh.selection.selected_strip_direction)); // flip selected strip direction
            cout << "Changed the selected_strip_direction to : " << np3dp->quad_mesh.selection.selected_strip_direction << endl;
            if (Strips::get_strip_from_fi(np3dp->quad_mesh.D, np3dp->quad_mesh.F, np3dp->quad_mesh.boundaryFaces, np3dp->quad_mesh.Emap, np3dp->quad_mesh.selection.selected_fi, np3dp->quad_mesh.selection.selected_strip_direction, np3dp->quad_mesh.selection.selected_strip, np3dp->quad_mesh.selection.selected_strip_vis, np3dp->quad_mesh.selection.strip_is_closed, vector<int>(), vector<int>())) {
                cout << "Traced strip with " << np3dp->quad_mesh.selection.selected_strip.sum() << " faces. Strip is closed : " << np3dp->quad_mesh.selection.strip_is_closed << endl;
                if (mark_helper_strips_viewer_index >=0) viewer.data_list[mark_helper_strips_viewer_index].clear();
            	display_selected_strip_data();
            }
            return true;
        }
        if (vertex_selection_is_on && np3dp->has_quad_mesh && !np3dp->quad_mesh.strip_networks.counterstrip_sequencies.empty()) { // change between different strip + counterstrip sequences.
			++np3dp->quad_mesh.strip_networks.selected_sequence_index;
            if (np3dp->quad_mesh.strip_networks.selected_sequence_index >= np3dp->quad_mesh.strip_networks.counterstrip_sequencies.size())
                np3dp->quad_mesh.strip_networks.selected_sequence_index = 0;

	        cout << "\nCurrent counterstrips sequence_index = 0 " << ", out of " << np3dp->quad_mesh.strip_networks.counterstrip_sequencies.size() - 1 << " total indices. Press SPACE to change selected sequence index" << endl;
	        cout << "The selected vertices are **" << np3dp->quad_mesh.strip_networks.strip_sequencies[np3dp->quad_mesh.strip_networks.selected_sequence_index].size() << "** strips apart in the blue direction." << endl;
            if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
            display_intermediate_strips_of_selected_vis(strip_sequences_viewer_index);
            return true;
        }
        return true;
    }
    return false;
}


void Gui::curl_reduction_one_step()
{
	std::cout << "--- Curl reduction step ---" << endl;
    std::cout << "Initial curl = " << vf->measure_U_curl(vf->matching) << endl;
    vf->project_U_curl_free(vf->matching);
    std::cout << "Final curl = " << vf->measure_U_curl(vf->matching) << endl;
    display_vector_field();
    update_mesh_colors();
}


void Gui::set_mesh_visual(const MatrixXd& V, const MatrixXi& F) {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    cout << "Displaying mesh, with nV = " << V.rows() << ", and nF = " << F.rows() << endl;

    // olga's color preferences :)
    viewer.data().set_colors(np3dp->mesh_base_color);
    viewer.core().background_color << 1., 1., 1., 1.;
    viewer.data().line_color << 173. / 255, 174. / 255, 103. / 255, 1.; //  np3dp->mesh_base_color[0], np3dp->mesh_base_color[1], np3dp->mesh_base_color[2], 0.0;// 
    viewer.data().line_width = 0.01; 

    // since mesh changes, remove any boundaries display
    if (mesh_boundaries_viewer_index >= 0) { viewer.data_list[mesh_boundaries_viewer_index].clear(); mesh_boundaries_viewer_index = -1; }
}


void Gui::update_mesh_colors() {
    show_vfConstraints_mesh_colors = meshColors == VfConstraints;
    show_face_degree_mesh_colors = meshColors == FaceDegree;
    show_quad_partition_map = meshColors == QuadPartitionMap;

    if (show_vfConstraints_mesh_colors && vf->constraint_fis.size() > 0 && !np3dp->has_quad_mesh) { display_fvf_constraints_face_colors(); }
    else if (show_face_degree_mesh_colors && np3dp->quad_mesh.V.size() > 0) { display_face_degree_mesh_colors(); }
    else if (show_quad_partition_map) { if (np3dp->Fmap.size() == np3dp->quad_mesh.F.rows()) display_per_face_data_on_quad_mesh(np3dp->quad_mesh.T, np3dp->quad_mesh.TF, np3dp->Fmap); }
    else { viewer.data().set_colors(np3dp->mesh_base_color); }// restore mesh colors 
}


void Gui::display_face_degree_mesh_colors() {
    int nF = np3dp->quad_mesh.T.rows();
    assert(np3dp->quad_mesh.TF.rows() == nF);
    MatrixXd face_colors(nF, 3);
    VectorXi degs(nF);
    for (int fi = 0; fi < np3dp->quad_mesh.T.rows(); ++fi) {
        int fi_original = np3dp->quad_mesh.TF[fi];
        degs[fi] = np3dp->quad_mesh.D[fi_original];
    }
    igl::parula(degs, true, face_colors);
    for (int fi = 0; fi < nF; ++fi)
        if (degs[fi] == 4)
            face_colors.row(fi) = np3dp->mesh_base_color;
    viewer.data().set_colors(face_colors);
}


void Gui::display_fvf_constraints_face_colors() {
    MatrixXd colors(np3dp->mesh.F.rows(), 3);
    vector<int> cfis(vf->constraint_fis.data(), vf->constraint_fis.data() + vf->constraint_fis.rows() * vf->constraint_fis.cols());
    for (int fi = 0; fi < np3dp->mesh.F.rows(); ++fi) {
        if (!Helpers::in_vector(cfis, fi)) colors.row(fi) = np3dp->mesh_base_color;
        else                               colors.row(fi) = RowVector3d(245.0 / 255.0, 171.0 / 255.0, 20 / 255.0);
    }
    viewer.data().set_colors(colors);
}


void Gui::prepare_viewer_index_data(int& index) {
    if (index < 0) {
        viewer.append_mesh();
        index = viewer.selected_data_index;
        viewer.data_list[index].show_custom_labels = true;
        viewer.selected_data_index = 0; // return to default
    }
    else {
        viewer.data_list[index].clear();
    }
}


void Gui::display_per_vertex_value(const VectorXd& values) {
    MatrixXd vertex_colors;
    igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, values, true, vertex_colors);
    viewer.data().set_colors(vertex_colors);
}


void Gui::display_Emap(const MatrixXd& Vquad, const MatrixXi& EVquad, const VectorXi& Emap)
{
    prepare_viewer_index_data(emap_viewer_index);
    viewer.data_list[emap_viewer_index].line_width = 3.0;

    vector<int> eis0, eis1, eis_boundary;
    for (int ei = 0; ei < Emap.size(); ++ei) {
        if (Emap[ei] == 0)      eis0.push_back(ei);
        else if (Emap[ei] == 1) eis1.push_back(ei);
        else                    eis_boundary.push_back(ei);
    }

    MatrixXd P0_start(eis0.size(), 3), P0_end(eis0.size(), 3);
    for (int i = 0; i < eis0.size(); ++i)
    {
        int ei = eis0[i];
        P0_start.row(i) = Vquad.row(EVquad(ei, 0));
        P0_end.row(i) = Vquad.row(EVquad(ei, 1));
    }
    MatrixXd P1_start(eis1.size(), 3), P1_end(eis1.size(), 3);
    for (int i = 0; i < eis1.size(); ++i)
    {
        int ei = eis1[i];
        P1_start.row(i) = Vquad.row(EVquad(ei, 0));
        P1_end.row(i) = Vquad.row(EVquad(ei, 1));
    }
    MatrixXd Pboundary_start(eis_boundary.size(), 3), Pboundary_end(eis_boundary.size(), 3);
    for (int i = 0; i < eis_boundary.size(); ++i)
    {
        int ei = eis_boundary[i];
        Pboundary_start.row(i) = Vquad.row(EVquad(ei, 0));
        Pboundary_end.row(i) = Vquad.row(EVquad(ei, 1));
    }

    viewer.data_list[emap_viewer_index].add_edges(P0_start, P0_end, RowVector3d(0, 0.5, 1));
    viewer.data_list[emap_viewer_index].add_edges(P1_start, P1_end, RowVector3d(0, 1, 0.5));
    viewer.data_list[emap_viewer_index].add_edges(Pboundary_start, Pboundary_end, RowVector3d(1, 0.5, 0));
}


void Gui::display_edge_maps()
{
    if (np3dp->has_quad_mesh && show_edge_maps)
    {
        if (show_ribs) {
            if (np3dp->piece_ID < 0 || np3dp->piece_ID >= np3dp->nP()) {
                if (np3dp->quad_mesh.ERibMap.rows() == np3dp->quad_mesh.Emap.rows())
                    display_Emap(np3dp->quad_mesh.V, np3dp->quad_mesh.EV, np3dp->quad_mesh.ERibMap);
            }
            else if (np3dp->has_partitioned_quad_mesh) {
                MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
                hedra::polygonal_edge_topology(np3dp->pieces[np3dp->piece_ID].D, np3dp->pieces[np3dp->piece_ID].F, EV, FE, EF, EFi, FEs, innerEdges);
                display_Emap(np3dp->pieces[np3dp->piece_ID].V, EV, np3dp->pieces[np3dp->piece_ID].ERibMap);
            }
            else if (emap_viewer_index >= 0) viewer.data_list[emap_viewer_index].clear();
        }
        else {
            if (np3dp->piece_ID == -1) {
                display_Emap(np3dp->quad_mesh.V, np3dp->quad_mesh.EV, np3dp->quad_mesh.Emap);
            }
            else if (np3dp->has_partitioned_quad_mesh) {
                MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
                hedra::polygonal_edge_topology(np3dp->pieces[np3dp->piece_ID].D, np3dp->pieces[np3dp->piece_ID].F, EV, FE, EF, EFi, FEs, innerEdges);
                display_Emap(np3dp->pieces[np3dp->piece_ID].V, EV, np3dp->pieces[np3dp->piece_ID].Emap);
            }
            else if (emap_viewer_index >= 0) viewer.data_list[emap_viewer_index].clear();
        }
    }
    else if (emap_viewer_index >= 0) viewer.data_list[emap_viewer_index].clear();
}


void Gui::mark_edges(const MatrixXd& V, const MatrixXi& EV, const vector<int>& eis, int& viewer_index, const RowVector3d& color, double line_width) {
    prepare_viewer_index_data(viewer_index);
    viewer.data_list[viewer_index].point_size = 8;
    viewer.data_list[viewer_index].line_width = line_width;

    MatrixXd P1(eis.size(), 3), P2(eis.size(), 3);
    for (int i = 0; i < eis.size(); ++i)
    {
        const int ei = eis[i];
        P1.row(i) = V.row(EV(ei, 0));
        P2.row(i) = V.row(EV(ei, 1));
    }
    viewer.data_list[viewer_index].add_edges(P1, P2, color);
}


void Gui::mark_edges(int& viewer_index, const RowVector3d& color, const MatrixXd& P1, const MatrixXd& P2, double line_width)
{
    prepare_viewer_index_data(viewer_index);
    viewer.data_list[viewer_index].point_size = 8;
    viewer.data_list[viewer_index].line_width = line_width;
    viewer.data_list[viewer_index].add_edges(P1, P2, color);
}


void Gui::mark_points(int& viewer_index, const RowVector3d& color, const MatrixXd& P) {
    viewer.data_list[viewer_index].point_size = 16;
    viewer.data_list[viewer_index].add_points(P, color);
}


/*void Gui::mark_fixed_vis()
{
    prepare_viewer_index_data(mark_fixed_pts_viewer_index);

    if (np3dp->FixedVis.rows() != np3dp->quad_mesh.V.rows()) { cerr << "V.rows()!= FixedVis.rows()  " << np3dp->quad_mesh.V.rows() << " != " << np3dp->FixedVis.rows() << endl; throw; }
    if (np3dp->FixedVis.sum() > 0)
    {
        MatrixXd P(np3dp->FixedVis.sum(), 3);
        int count = 0;
        for (int vi = 0; vi < np3dp->quad_mesh.V.rows(); ++vi) {
            if (np3dp->FixedVis[vi]) {
                P.row(count) = np3dp->quad_mesh.V.row(vi);
                ++count;
            }
        }
        mark_points(mark_fixed_pts_viewer_index, RowVector3d(0.55, 0.05, 0.85), P);
    }
}*/


void Gui::mark_quad_edges() {
    prepare_viewer_index_data(quad_edges_viewer_index);
    viewer.data_list[quad_edges_viewer_index].line_width = 8.0;
    const MatrixXd& Vquad = np3dp->quad_mesh.V;
    const MatrixXi& EVquad = np3dp->quad_mesh.EV;
    const RowVector3d color(0.0, 0.5, 1.0);

    MatrixXd P1(EVquad.rows(), 3), P2(EVquad.rows(), 3);
    for (int ei = 0; ei < EVquad.rows(); ++ei) {
        P1.row(ei) = Vquad.row(EVquad(ei, 0));
        P2.row(ei) = Vquad.row(EVquad(ei, 1));
    }
    viewer.data_list[quad_edges_viewer_index].add_edges(P1, P2, color);
}


void Gui::display_mesh_boundaries(const MatrixXd& V, const MatrixXi& F) {
    const RowVector3d color(0.9, 0, 0.5);

    if (!np3dp->has_quad_mesh)
    {
        MatrixXi EV, FE, EF;
        igl::edge_topology(V, F, EV, FE, EF);

        vector<int> boundary_eis;
        for (int ei = 0; ei < EF.rows(); ++ei)
            if (EF(ei, 0) < 0 || EF(ei, 1) < 0)
                boundary_eis.push_back(ei);

        mark_edges(V, EV, boundary_eis, mesh_boundaries_viewer_index, color);
    }
    else
    {
        // for now just show boundary vertices
        prepare_viewer_index_data(mesh_boundaries_viewer_index);
        viewer.data_list[mesh_boundaries_viewer_index].point_size = 8;
        for (int vi = 0; vi < np3dp->quad_mesh.boundaryVertices.rows(); ++vi)
            if (np3dp->quad_mesh.boundaryVertices[vi] == 1)
                viewer.data_list[mesh_boundaries_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), color);
    }
}


void Gui::display_vf_singularities() { // VF singularities are on the mesh VERTICES
    prepare_viewer_index_data(vf_sings_viewer_index);
    viewer.data_list[vf_sings_viewer_index].point_size = 55;

    if (vf->U.size() == 2 * np3dp->mesh.F.rows()) {
        MatrixXd cartesian_field, raw_field;
		MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesian_field, vf->U.head(np3dp->mesh.F.rows()));
		directional::representative_to_raw(np3dp->mesh.faceNormals, cartesian_field, vf->fieldDegree, raw_field);
  
		directional::CartesianField field(vf->ftb);
		field.N = vf->fieldDegree;
		field.fieldType = directional::fieldTypeEnum::RAW_FIELD;
		field.set_extrinsic_field(raw_field);
		directional::principal_matching(field);

        VectorXi sing_vertices = field.singLocalCycles; // remember, this also contains all boundary vertices!
    	VectorXi sing_indices = field.singIndices;

        RowVector3d color(1, 0, 0);
        for (int i = 0; i < sing_vertices.size(); ++i) {
            int vi = sing_vertices[i];
            viewer.data_list[vf_sings_viewer_index].add_points(np3dp->mesh.V.row(vi), color);
        }
        cout << endl;
    }
}


void Gui::display_quad_mesh_vertex_degrees() {
    prepare_viewer_index_data(quad_vertex_degrees_viewer_index);
    viewer.data_list[quad_vertex_degrees_viewer_index].point_size = 15;
    MatrixXd colors;
    igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, np3dp->quad_mesh.vertexDegree, true, colors);
    int count = 0;
    for (int vi = 0; vi < np3dp->quad_mesh.V.rows(); ++vi) {
        if (np3dp->quad_mesh.vertexDegree[vi] != 4 && np3dp->quad_mesh.boundaryVertices[vi] != 1){
	        viewer.data_list[quad_vertex_degrees_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), colors.row(vi));
        	++count;
        }
    }
    cout << "Number of irregular vertices = " << count << endl;
}


void Gui::display_vector_field()
{
    if (show_vector_field)
    {
        prepare_viewer_index_data(vf_viewer_index);
        prepare_viewer_index_data(vf_rotated_viewer_index);
        int nF = np3dp->mesh.F.rows();
        if (vf->U.rows() == 2 * np3dp->mesh.F.rows()) {
        	vf->display_vector_field(viewer, vf_viewer_index, vf->U.head(nF), 0);
            vf->save_vf_to_txt("vf.txt", vf->U.head(nF));
            vf->display_vector_field(viewer, vf_rotated_viewer_index, vf->U.tail(nF), 1);
            vf->save_vf_to_txt("vf_rotated.txt", vf->U.tail(nF));
        }
    }
    if (show_vf_sings) display_vf_singularities();
}


void Gui::display_paths() {
    if (paths_viewer_index >= 0) viewer.data_list[paths_viewer_index].clear();

    if (np3dp->piece_ID >= 0 && np3dp->piece_ID < np3dp->nP()) { // only display ID piece
        np3dp->pieces[np3dp->piece_ID].pathsCollection.display(viewer, paths_viewer_index, true);
    }
    else { // display all pieces
        for (int id = 0; id < np3dp->nP(); ++id)
            np3dp->pieces[id].pathsCollection.display(viewer, paths_viewer_index, false);
    }
}


void Gui::display_UV_coords(const MatrixXd& VMeshCut, const MatrixXi& FMeshCut, const MatrixXd& UVcoords)
{
    set_mesh_visual(VMeshCut, FMeshCut);
    viewer.data().set_uv(UVcoords);

    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>  R, G, B;
    Helpers::setup_line_texture(R, G, B, np3dp->texture_scale, np3dp->texture_line_width);
    viewer.data().set_texture(R, G, B);
    viewer.data().show_texture = true;
}


void Gui::display_partitioned_mesh() {
    meshColors = QuadPartitionMap;
    display_per_face_data_on_quad_mesh(np3dp->quad_mesh.T, np3dp->quad_mesh.TF, np3dp->Fmap);
    vector<int> eis;
    for (int ei = 0; ei < np3dp->PartitioningCuts.rows(); ++ei) if (np3dp->PartitioningCuts[ei] == 1) eis.push_back(ei);
    mark_edges(np3dp->quad_mesh.V, np3dp->quad_mesh.EV, eis, cuts_on_edges_viewer_index, RowVector3d(1, 0, 0));
    display_edge_maps();
}


void Gui::draw_partitioning_cuts() {
    VectorXi currentPartitioningCuts;
    Partitioning::propagate_cuts_on_quad_mesh(np3dp.get(), currentPartitioningCuts);
    vector<int> eis;
    for (int ei = 0; ei < currentPartitioningCuts.rows(); ++ei) if (currentPartitioningCuts[ei] == 1) eis.push_back(ei);
    mark_edges(np3dp->quad_mesh.V, np3dp->quad_mesh.EV, eis, cuts_on_edges_viewer_index, RowVector3d(1, 0, 0), 8.0);

    Helpers::write_matrix_to_txt(np3dp->quad_mesh.get_lines_matrix_from_cuts(currentPartitioningCuts), DATA_PATH + np3dp->output_folder + "current_cuts.txt");

}


void Gui::display_selected_strip_data()
{
    // face colors
    int nF = np3dp->quad_mesh.T.rows();
    assert(np3dp->quad_mesh.TF.rows() == nF);
    MatrixXd face_colors = np3dp->mesh_base_color.replicate(nF, 1);

    const RowVector3d strip_color = np3dp->quad_mesh.selection.selected_strip_direction == 0 ? np3dp->dominant_strip_color : np3dp->subdominant_strip_color;

    for (int fi=0; fi<np3dp->quad_mesh.T.rows(); ++fi)
        if (np3dp->quad_mesh.selection.selected_strip[np3dp->quad_mesh.TF[fi]])
            face_colors.row(fi) = strip_color;
    viewer.data().set_colors(face_colors);

    // points
    {
        prepare_viewer_index_data(mark_selected_viewer_index);
        viewer.data_list[mark_selected_viewer_index].point_size = 10;
        for (int i = 0; i < 2; ++i) { // once for each side of the strip
            const vector<int>& vis = np3dp->quad_mesh.selection.selected_strip_vis[i];
            RowVector3d color = i == 0 ? RowVector3d(1.0, 0.5, 0.0) : RowVector3d(0.0, 1.0, 0.5);
            int count = 0;
            for (int vi : vis) {
                viewer.data_list[mark_selected_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), color);
                //viewer.data_list[mark_selected_viewer_index].add_label(np3dp->quad_mesh.V.row(vi), std::to_string(count));
                ++count;
            }
        }
    }

    //// arrows
    //{
    //    viewer.data_list[mark_selected_viewer_index].line_width = 15;
    //    const vector<int>& vis_A = np3dp->strip_vis[0];
    //    RowVector3d color(0.32, 0.1, 0.88);
    //    for (int k = 0; k < vis_A.size() - 1; ++k) {
    //        const int vi1 = vis_A[k];
    //        const RowVector3d& v1 = np3dp->quad_mesh.V.row(vi1);
    //        const int vi2 = vis_A[k + 1];
    //        const RowVector3d& v2 = np3dp->quad_mesh.V.row(vi2);
    //        RowVector3d ei_vec = v1 - v2;
    //        //viewer.data_list[mark_selected_viewer_index].add_edges(v1 + ei_vec * 0.35, v2 - ei_vec * 0.1, color);
    //        viewer.data_list[mark_selected_viewer_index].add_points(v2 - ei_vec * 0.2, color);
    //    }
    //}

    // centerlines
    {
        //display_strips_centerlines(mark_helper_strips_viewer_index, np3dp->strip, strip_color * 0.5, np3dp->selected_strip_direction);
    }
}


void Gui::display_strips_centerlines(int& viewer_index, const VectorXi& strip, const RowVector3d& color, int direction)
{
    if (viewer_index < 0) prepare_viewer_index_data(viewer_index);
    viewer.data_list[viewer_index].point_size = 10;
    viewer.data_list[viewer_index].line_width = 5;

    int other_direction = static_cast<int>(!static_cast<bool>(direction));

    auto get_fcen = [](const MatrixXi& F, const MatrixXd& V, int fi, const VectorXi& D)->RowVector3d {
        RowVector3d cen(0.0, 0.0, 0.0);
        for (int k = 0; k < D[fi]; ++k)
            cen += V.row(F(fi, k));
        cen = cen.cwiseProduct(RowVector3d(1.0 / D[fi], 1.0 / D[fi], 1.0 / D[fi]));
        return cen;
    };

    const MatrixXi& F = np3dp->quad_mesh.F;
    const MatrixXd& V = np3dp->quad_mesh.V;
    VectorXi passed; passed.setZero(F.rows());

    vector<RowVector3d> pts;
    vector<RowVector3d> P1s, P2s;

    for (int fi = 0; fi < F.rows(); ++fi) {
        if (strip[fi]) {
            if (!passed[fi]) {
                passed[fi] = 1;

                RowVector3d cen = get_fcen(F, V, fi, np3dp->quad_mesh.D);
                //cen += N.row(fi) * 0.01;

                // add point at center of face
                viewer.data_list[viewer_index].add_points(cen, color);

                // find neighbor that has not been passed, and add line to that neighbor
                for (int k = 0; k < np3dp->quad_mesh.D[fi]; ++k)
                {
                    int ei = np3dp->quad_mesh.FE(fi, k);
                    if (np3dp->quad_mesh.Emap[ei] == other_direction)
                    {
                        int nf = np3dp->quad_mesh.EF(ei, 0) == fi ? np3dp->quad_mesh.EF(ei, 1) : np3dp->quad_mesh.EF(ei, 0);
                        if (nf >= 0) {
                            if (strip[nf]) {
                                if (!passed[nf]) {
                                    RowVector3d other_cen = get_fcen(F, V, nf, np3dp->quad_mesh.D);
                                    //other_cen += N.row(nf) * 0.01;

                                    // add point at center of face
                                    pts.push_back(other_cen);

                                    // add line connecting two faces
                                    P1s.push_back(cen);
                                    P2s.push_back(other_cen);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    MatrixXd Points(pts.size(), 3);
    for (int i = 0; i < pts.size(); ++i)
        Points.row(i) = pts[i];

    MatrixXd P1(P1s.size(), 3), P2(P2s.size(), 3);
    for (int i = 0; i < P1s.size(); ++i)
    {
        P1.row(i) = P1s[i];
        P2.row(i) = P2s[i];
    }

    viewer.data_list[viewer_index].add_points(Points, color);

    // add line connecting two faces
    viewer.data_list[viewer_index].add_edges(P1, P2, color);
}


void Gui::display_markings(const vector<Vector3d>& positions)
{
    prepare_viewer_index_data(position_markers_viewer_index);

    MatrixXd markings(positions.size(), 3);
    for (int i = 0; i < positions.size(); ++i)
        markings.row(i) = RowVector3d(positions[i]);
    viewer.data_list[position_markers_viewer_index].point_size = 15.0;
    viewer.data_list[position_markers_viewer_index].add_points(markings, RowVector3d(0.9, 0.05, 0.25));
}


void Gui::display_selected_vis(const RowVector3d& color)
{
    prepare_viewer_index_data(mark_selected_vertices_viewer_index);
    for (int vi : np3dp->quad_mesh.selection.selected_vis)
        viewer.data_list[mark_selected_vertices_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), color);
}


void Gui::display_singularity_alignments(int& viewer_index, const MatrixXd& P1, const MatrixXd& P2, const vector<int>& vis, const RowVector3d& color, const RowVector3d& other_color)
{
    prepare_viewer_index_data(viewer_index);

    // show path to singularities
    viewer.data_list[viewer_index].point_size = 25;
    viewer.data_list[viewer_index].line_width = 20.0;
    viewer.data_list[viewer_index].add_edges(P1, P2, color);

    // show aligned singularities
    MatrixXd P; P.setZero(vis.size(), 3);
    for (int i = 0; i<vis.size(); ++i)
        P.row(i) = np3dp->quad_mesh.V.row(vis[i]);
	viewer.data_list[viewer_index].add_points(P, color);

    // also show tracked vis that have not been added to some alignment
    vector<int> other_vis;
    for (pair<int, int>& pair : np3dp->quad_mesh.vis_to_check_for_alignments){
        int v1 = pair.first;
        int v2 = pair.second;
        if (!Helpers::in_vector(vis, v1)) other_vis.push_back(v1);
        if (!Helpers::in_vector(vis, v2)) other_vis.push_back(v2);
    }

    P.setZero(other_vis.size(), 3);
    for (int i = 0; i<other_vis.size(); ++i)
        P.row(i) = np3dp->quad_mesh.V.row(other_vis[i]);
	viewer.data_list[viewer_index].add_points(P, other_color);
}


void Gui::display_blocked_strips()
{
	const StripNetworks& sn = np3dp->quad_mesh.strip_networks;

    int nF = np3dp->quad_mesh.T.rows();
    MatrixXd face_colors = np3dp->mesh_base_color.replicate(nF, 1);

    for (int fi = 0; fi < np3dp->quad_mesh.T.rows(); ++fi) {
        int quad_fi = np3dp->quad_mesh.TF[fi];
        for (int si = 0; si < sn.blockedStrips.size(); ++si) {
            if (sn.blockedStrips[si]) {
                if (np3dp->quad_mesh.strip_networks.strip_index_to_properties(si).strip[quad_fi]) { // if this face is part of a blocked string
                    face_colors.row(fi) = np3dp->mesh_base_color * 0.6;
                    break;
                }
            }

        }
    }
    viewer.data().set_colors(face_colors);
}


void Gui::display_intermediate_strips_of_selected_vis(int& viewer_index)
{
    // --- get display strips
	const StripNetworks& sn = np3dp->quad_mesh.strip_networks;

    MatrixXi strips;
    if (sn.strip_sequencies.size() > sn.selected_sequence_index )
    {
    	strips.setZero(sn.strip_sequencies[sn.selected_sequence_index].size(), np3dp->quad_mesh.F.rows());
	    for (int i = 0; i < sn.strip_sequencies[sn.selected_sequence_index].size(); ++i) {
	        int si = sn.strip_sequencies[sn.selected_sequence_index][i];
	        strips.row(i) = RowVectorXi(sn.Dgraph[sn.StoG[si]].strip);
	    }
    }

    MatrixXi counterstrips;
    counterstrips.setZero(sn.counterstrip_sequencies[sn.selected_sequence_index].size(), np3dp->quad_mesh.F.rows());
    for (int i = 0; i < sn.counterstrip_sequencies[sn.selected_sequence_index].size(); ++i) {
        int si = sn.counterstrip_sequencies[sn.selected_sequence_index][i];
        counterstrips.row(i) = RowVectorXi(sn.cDgraph[sn.StoG[si]].strip);
    }

    MatrixXi blocked_strips;
    blocked_strips.setZero(sn.blockedStrips.sum(), np3dp->quad_mesh.F.rows());
    int count = 0;
    for (int si = 0; si <sn.blockedStrips.size(); ++si) {
        if (sn.blockedStrips[si] == 1){
            const StripGraph& G = sn.StoD[si] == 0 ? sn.Dgraph : sn.cDgraph;
            blocked_strips.row(count) = RowVectorXi(G[sn.StoG[si]].strip);
            ++count;
        }
    }

    // --- create face colors
    int nF = np3dp->quad_mesh.T.rows();
    MatrixXd face_colors = np3dp->mesh_base_color.replicate(nF, 1);
    for (int fi = 0; fi < nF; ++fi) {
        if (strips.rows() > 0 && counterstrips.rows() > 0){ 
            if (strips.col(np3dp->quad_mesh.TF[fi]).maxCoeff() == 1 && counterstrips.col(np3dp->quad_mesh.TF[fi]).maxCoeff() == 1) { // face in both strips and counterstrips
                face_colors.row(fi) = 0.5 * (np3dp->dominant_strip_color + np3dp->subdominant_strip_color);
                goto ctn;
            }
        }
        if (strips.rows() > 0) {
	        if (strips.col(np3dp->quad_mesh.TF[fi]).maxCoeff() == 1){ // face in strips
                face_colors.row(fi) = np3dp->dominant_strip_color;
                continue; // dominant strips should not be made darker if they are blocked (because their collapse would not change alignment)
	        }
        }
        if (counterstrips.rows() > 0) {
			if (counterstrips.col(np3dp->quad_mesh.TF[fi]).maxCoeff() == 1){ // face in counterstrips
                face_colors.row(fi) = np3dp->subdominant_strip_color;
			}
        }
        ctn:;
        if (sn.blockedStrips.sum() > 0){
	        if (blocked_strips.col(np3dp->quad_mesh.TF[fi]).maxCoeff() == 1){
                face_colors.row(fi) *= 0.3;
	        }
        }
    }
    viewer.data().set_colors(face_colors);


    // --- display vis sequences
    prepare_viewer_index_data(viewer_index);

    const vector<int>& vis = sn.vis_sequences[sn.selected_sequence_index];
    MatrixXd P1(vis.size()-1, 3), P2(vis.size()-1, 3);
    for (int i=1; i<vis.size(); ++i){
	    P1.row(i-1) = np3dp->quad_mesh.V.row(vis[i-1]);
        P2.row(i-1) = np3dp->quad_mesh.V.row(vis[i]);
    }
    RowVector3d color(0.84, 0.48, 0.32);
    viewer.data_list[viewer_index].add_edges(P1, P2, color);
    // P1.conservativeResize(P1.rows()+1, 3);
    // P1.row(P1.rows()-1) = np3dp->quad_mesh.V.row(vis.back());
    // viewer.data_list[viewer_index].add_points(P1, color * 0.75);

	viewer.data_list[viewer_index].label_size = 2;
	viewer.data_list[viewer_index].clear_labels();
    for (int k=0; k<vis.size(); ++k)
    {
		int vi = vis[k];
        viewer.data_list[mark_selected_vertices_viewer_index].add_label(np3dp->quad_mesh.V.row(vi), std::to_string(k));
    }
}


void Gui::update_visualization_after_topology_edits()
{
	if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	if (show_quad_vertex_degrees) display_quad_mesh_vertex_degrees();
    if (cuts_on_edges_viewer_index >= 0) draw_partitioning_cuts();
    if (mark_fixed_pts_viewer_index >=0) viewer.data_list[mark_fixed_pts_viewer_index].clear();
    if (strip_sequences_viewer_index >=0) viewer.data_list[strip_sequences_viewer_index].clear();

    set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
    display_edge_maps();
    display_selected_vis();
    update_mesh_colors();
}


void Gui::invalidate(bool keep_vf_params) {
	np3dp = make_shared<Np3dpContext>(np3dp->data_folder);

    if (!keep_vf_params)
		vf = make_shared<VectorField>(np3dp);
    else
        vf = make_shared<VectorField>(np3dp, vf->time_steps, vf->vfGenerationType,  vf->f_fraction, vf->smoothness_coeff, vf->alignment_coeff, vf->orthogonality_coeff, vf->unit_coeff, vf->add_boundary_constraints, vf->curvature_absolute_values);

    meshColors = VfConstraints;
    show_vector_field = true;
    show_vf_sings = true;

    set_mesh_visual(np3dp->mesh.V, np3dp->mesh.F);
    display_vector_field();
    update_mesh_colors();

    np3dp->piece_ID = -1;

    if (vf_viewer_index >= 0) { viewer.data_list[vf_viewer_index].clear(); vf_viewer_index = -1; }
    if (vf_rotated_viewer_index >= 0) { viewer.data_list[vf_rotated_viewer_index].clear(); vf_viewer_index = -1; }
    if (vf_sings_viewer_index >= 0) { viewer.data_list[vf_sings_viewer_index].clear(); vf_sings_viewer_index = -1; }
    if (vf_streamlines_viewer_index >= 0) { viewer.data_list[vf_streamlines_viewer_index].clear(); vf_streamlines_viewer_index = -1; }
	if (mesh_boundaries_viewer_index >= 0) { viewer.data_list[mesh_boundaries_viewer_index].clear(); mesh_boundaries_viewer_index = -1; }
    if (paths_viewer_index >= 0) { viewer.data_list[paths_viewer_index].clear(); paths_viewer_index = -1; }
    if (strip_sequences_viewer_index >= 0) { viewer.data_list[strip_sequences_viewer_index].clear(); strip_sequences_viewer_index = -1; }
    if (quad_edges_viewer_index >= 0) { viewer.data_list[quad_edges_viewer_index].clear(); quad_edges_viewer_index = -1; }
    if (quad_vertex_degrees_viewer_index >= 0) { viewer.data_list[quad_vertex_degrees_viewer_index].clear(); quad_vertex_degrees_viewer_index = -1; }
    if (emap_viewer_index >= 0) { viewer.data_list[emap_viewer_index].clear(); emap_viewer_index = -1; }
    if (cuts_on_edges_viewer_index >= 0) { viewer.data_list[cuts_on_edges_viewer_index].clear(); cuts_on_edges_viewer_index = -1; }
    if (mark_selected_viewer_index >= 0) { viewer.data_list[mark_selected_viewer_index].clear(); mark_selected_viewer_index = -1; }
    if (mark_selected_vertices_viewer_index >= 0) { viewer.data_list[mark_selected_vertices_viewer_index].clear(); mark_selected_vertices_viewer_index = -1; }
    if (mark_fixed_pts_viewer_index >= 0) { viewer.data_list[mark_fixed_pts_viewer_index].clear(); mark_fixed_pts_viewer_index = -1; }
    if (mark_helper_strips_viewer_index >= 0) { viewer.data_list[mark_helper_strips_viewer_index].clear(); mark_helper_strips_viewer_index = -1; }
    if (alignment_of_selected_vertices_viewer_index >= 0) { viewer.data_list[alignment_of_selected_vertices_viewer_index].clear(); alignment_of_selected_vertices_viewer_index = -1; }
    if (position_markers_viewer_index >= 0) { viewer.data_list[position_markers_viewer_index].clear(); position_markers_viewer_index = -1; }
    if (piece_labels_viewer_index >= 0) { viewer.data_list[piece_labels_viewer_index].clear(); piece_labels_viewer_index = -1; }

	cout << "Invalidated all" << endl << endl;
}


void Gui::serialize() const
{
    cout << endl;

    // empty serialize folder
    Helpers::remove_directory_contents(DATA_PATH + np3dp->serialize_folder);

    np3dp->serialize();
    vf->serialize();


    igl::serialize(show_mesh_boundaries, DATA_PATH + np3dp->serialize_folder + "gui.show_mesh_boundaries.igl");
    igl::serialize(show_vf_sings, DATA_PATH + np3dp->serialize_folder + "gui.show_vf_sings.igl");
    igl::serialize(show_vector_field, DATA_PATH + np3dp->serialize_folder + "gui.show_vector_field.igl");
    igl::serialize(show_paths, DATA_PATH + np3dp->serialize_folder + "gui.show_paths.igl");
    igl::serialize(show_paths_rotated, DATA_PATH + np3dp->serialize_folder + "gui.show_paths_rotated.igl");
    igl::serialize(show_quad_vertex_degrees, DATA_PATH + np3dp->serialize_folder + "gui.show_quad_vertex_degrees.igl");
    igl::serialize(show_helper_strips, DATA_PATH + np3dp->serialize_folder + "gui.show_helper_strips.igl");
    igl::serialize(show_helper_counter_strips, DATA_PATH + np3dp->serialize_folder + "gui.show_helper_counter_strips.igl");
    igl::serialize(show_edge_maps, DATA_PATH + np3dp->serialize_folder + "gui.show_edge_maps.igl");
    igl::serialize(show_ribs, DATA_PATH + np3dp->serialize_folder + "gui.show_ribs.igl");
    igl::serialize(show_vfConstraints_mesh_colors, DATA_PATH + np3dp->serialize_folder + "gui.show_vfConstraints_mesh_colors.igl");
    igl::serialize(show_face_degree_mesh_colors, DATA_PATH + np3dp->serialize_folder + "gui.show_face_degree_mesh_colors.igl");
    igl::serialize(show_quad_partition_map, DATA_PATH + np3dp->serialize_folder + "gui.show_quad_partition_map.igl");
	igl::serialize(show_singularities_alignment, DATA_PATH + np3dp->serialize_folder + "gui.show_singularities_alignment.igl");
	
    std::cout << "Serialized all data" << endl << endl;
}


void Gui::deserialize()
{
    invalidate();
    np3dp->deserialize();

    bool do_vf_setup = !np3dp->has_quad_mesh && !np3dp->has_uv_coords;
    vf->deserialize(do_vf_setup);

    igl::deserialize(show_mesh_boundaries, DATA_PATH + np3dp->serialize_folder + "gui.show_mesh_boundaries.igl");
    igl::deserialize(show_vf_sings, DATA_PATH + np3dp->serialize_folder + "gui.show_vf_sings.igl");
    igl::deserialize(show_vector_field, DATA_PATH + np3dp->serialize_folder + "gui.show_vector_field.igl");
    igl::deserialize(show_paths, DATA_PATH + np3dp->serialize_folder + "gui.show_paths.igl");
    igl::deserialize(show_paths_rotated, DATA_PATH + np3dp->serialize_folder + "gui.show_paths_rotated.igl");
    igl::deserialize(show_quad_vertex_degrees, DATA_PATH + np3dp->serialize_folder + "gui.show_quad_vertex_degrees.igl");
    igl::deserialize(show_helper_strips, DATA_PATH + np3dp->serialize_folder + "gui.show_helper_strips.igl");
    igl::deserialize(show_helper_counter_strips, DATA_PATH + np3dp->serialize_folder + "gui.show_helper_counter_strips.igl");
    igl::deserialize(show_edge_maps, DATA_PATH + np3dp->serialize_folder + "gui.show_edge_maps.igl");
    igl::deserialize(show_ribs, DATA_PATH + np3dp->serialize_folder + "gui.show_ribs.igl");
    igl::deserialize(show_vfConstraints_mesh_colors, DATA_PATH + np3dp->serialize_folder + "gui.show_vfConstraints_mesh_colors.igl");
    igl::deserialize(show_face_degree_mesh_colors, DATA_PATH + np3dp->serialize_folder + "gui.show_face_degree_mesh_colors.igl");
    igl::deserialize(show_quad_partition_map, DATA_PATH + np3dp->serialize_folder + "gui.show_quad_partition_map.igl");
	igl::deserialize(show_singularities_alignment, DATA_PATH + np3dp->serialize_folder + "gui.show_singularities_alignment.igl");

	// --- Display deserialized data
	np3dp->has_uv_coords = false; // integration should be setup again
	np3dp->has_paths_on_pieces = false; // paths are not serialized

    if (!np3dp->has_quad_mesh) {
        display_vector_field();
        update_mesh_colors();
    }
    else if (np3dp->has_quad_mesh && !np3dp->has_partitioned_quad_mesh)
    {
        set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
        display_quad_mesh_vertex_degrees();
        display_edge_maps();
        // display_UV_coords(np3dp->quad_mesh.V, np3dp->quad_mesh.T, np3dp->quad_mesh.UVcoords);
        meshColors = FaceDegree;
        update_mesh_colors();
        np3dp->quad_mesh.save_data();
    }
    else if (np3dp->has_quad_mesh && np3dp->has_partitioned_quad_mesh)
    {
        set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
        display_edge_maps();
        display_partitioned_mesh();
        np3dp->quad_mesh.save_data();
    }

    std::cout << "Deserialized all data" << endl << endl;
}


void Gui::draw_menu_boxes()
{

	menu.callback_draw_viewer_menu = [&]() {
	menu.draw_viewer_menu();

	int first_win_size = 237;

	// --- Utilities window
	{
	    ImGui::SetNextWindowPos(ImVec2(0, 762));
	    ImGui::SetNextWindowSize(ImVec2(first_win_size, 0));
	    ImGui::Begin("Utilities");
	    if (ImGui::Button("(Q) Deserialize", ImVec2(115, 0))) deserialize();
	    ImGui::SameLine(126);
	    if (ImGui::Button("Serialize", ImVec2(-1, 0))) serialize();
	    if (ImGui::Button("Invalidate all", ImVec2(-1, 0))) invalidate(true);

	    // Display options
	    {
	        ImGui::Text(""); ImGui::Text("--- Display Options ---");
	        if (ImGui::Combo("Mesh colors", (int*)(&meshColors), "None\0VfConstraints\0FaceDegree\0QuadPartitionMap\0"))
	            update_mesh_colors();
	        if (ImGui::Checkbox("VF", &show_vector_field)) {
	            if (show_vector_field && vf->is_setup) display_vector_field();
	            else {
	                if (vf_viewer_index >= 0) { viewer.data_list[vf_viewer_index].clear(); viewer.data_list[vf_rotated_viewer_index].clear(); }
	            }
	        }
	        ImGui::SameLine(65);
	        if (ImGui::Checkbox("Sings", &show_vf_sings)) {
	            if (show_vf_sings && vf->is_setup) display_vf_singularities();
	            else if (vf_sings_viewer_index >= 0) { viewer.data_list[vf_sings_viewer_index].clear(); }
	        }
	        if (ImGui::InputDouble("vf_sizeRatio", &vf->sizeRatio)) if (show_vector_field && vf->is_setup) display_vector_field();
	        if (ImGui::InputInt("vf_sparsity", &vf->sparsity)) if (show_vector_field && vf->is_setup) display_vector_field();

	        if (ImGui::Checkbox("Show Mesh Boundaries", &show_mesh_boundaries)) {
	            if (show_mesh_boundaries) display_mesh_boundaries(np3dp->mesh.V, np3dp->mesh.F);
	            else if (mesh_boundaries_viewer_index >= 0) viewer.data_list[mesh_boundaries_viewer_index].clear();
	        }
	        if (ImGui::Checkbox("Show Paths", &show_paths)) {
	            display_paths();
	        }
	        if (np3dp->has_quad_mesh)
	        {
	            if (ImGui::Checkbox("Show Quad mesh vertex degrees", &show_quad_vertex_degrees)) {
	                if (show_quad_vertex_degrees && np3dp->quad_mesh.vertexDegree.size() == np3dp->quad_mesh.V.rows()) display_quad_mesh_vertex_degrees();
	                else if (quad_vertex_degrees_viewer_index >= 0) viewer.data_list[quad_vertex_degrees_viewer_index].clear();
	            }
	        }


	        if (ImGui::Checkbox("Show edges", &show_edge_maps))
	            display_edge_maps();
	        ImGui::SameLine(140);
	        if (ImGui::Checkbox("Show ribs", &show_ribs))
	            display_edge_maps();


	        if (ImGui::Checkbox("Show strips", &show_helper_strips))
	            display_intermediate_strips_of_selected_vis(strip_sequences_viewer_index);
	        ImGui::SameLine(140);
	        if (ImGui::Checkbox("Show counter-strips", &show_helper_strips))
	            display_intermediate_strips_of_selected_vis(strip_sequences_viewer_index);
	        
	        if (ImGui::Checkbox("Show singularities alignments", &show_singularities_alignment)){
	            MatrixXd P1, P2; vector<int> vis;
	            for (const auto& a : np3dp->quad_mesh.aligned_vertices)
	                a.append_visualization_data(np3dp->quad_mesh.V, P1, P2, vis);
	            display_singularity_alignments(strip_sequences_viewer_index, P1, P2, vis);
	        }
	    }
	    ImGui::End();
	}


	ImGui::SetNextWindowPos(ImVec2(first_win_size, 0));
	ImGui::SetNextWindowSize(ImVec2(420, 0));
	ImGui::Begin("Np-3DP");

	// --- Vector field generation
	{
	    ImGui::Text("--- Vector field ---");

	    // parameters
	    if (ImGui::Combo("Align with", (int*)(&vf->vfGenerationType), "Curvature\0Boundaries\0ConstraintsFromFile\0Random\0"))
	    {
	        vf->vfGenerationType = vf->vfGenerationType;
	        if (vf->vfGenerationType == VFGenerationType::CurvatureAligned) {
	            vf->alignment_coeff = 100;
	            vf->f_fraction = 0.2;
	        }
	        else if (vf->vfGenerationType == VFGenerationType::Random) {
	            vf->f_fraction = 0.02;
	        }
	        else { // boundaries algined / constraints from file 
	        }

	        invalidate(true);
	    }
	    if (ImGui::Combo("Time steps type", (int*)(&vf->time_steps), "Explicit\0Implicit\0")) invalidate(true);

	    // if (vf->vfGenerationType == VFGenerationType::CurvatureAligned)
	    //     if (ImGui::InputInt("Curvature confidence type: 0/1/2", &np3dp->curvature_confidence_type))
	    //         invalidate(true);
	    if (vf->vfGenerationType == VFGenerationType::CurvatureAligned){
	        if (ImGui::Checkbox("Consider curvature absolute values", &vf->curvature_absolute_values)){
		        invalidate(true);
	        }
	    }


	    ImGui::Text("Set following constraints in both directions");
	    ImGui::Checkbox("Alignment", &vf->set_alignment_constraints_in_both_directions);
	    ImGui::SameLine(140);
	    ImGui::Checkbox("Smoothness", &vf->set_smoothness_constraint_in_both_directions);
	    ImGui::SameLine(280);
	    ImGui::Checkbox("Unit", &vf->set_unit_constraint_in_both_directions);


	    if (ImGui::InputDouble("Smoothness coefficient", &vf->smoothness_coeff)) { if (vf->is_setup && vf->time_steps == TimeSteps::Implicit) vf->prepare_implicit_steps_solver(); }
	    if (ImGui::InputDouble("Alignment coefficient", &vf->alignment_coeff)) { if (vf->is_setup && vf->time_steps == TimeSteps::Implicit) vf->prepare_implicit_steps_solver(); }
	    if (ImGui::InputDouble("Orthogonality coefficient", &vf->orthogonality_coeff)) { if (vf->is_setup && vf->time_steps == TimeSteps::Implicit) vf->prepare_implicit_steps_solver(); }
	    if (ImGui::InputDouble("Unitize coefficient", &vf->unit_coeff)) { if (vf->is_setup && vf->time_steps == TimeSteps::Implicit) vf->prepare_implicit_steps_solver(); }
	    ImGui::InputDouble("Coefficients multiplier per iteration", &vf->decrease_mult);


	    if (vf->vfGenerationType == VFGenerationType::CurvatureAligned || vf->vfGenerationType == VFGenerationType::Random) { if (ImGui::InputDouble("% faces", &vf->f_fraction)) invalidate(true); }
	    
	    // if (ImGui::Checkbox("Add boundary constraints", &vf->add_boundary_constraints)) invalidate(true);
	    ImGui::InputInt("Max iterations", &vf->max_iterations);

	    /*if (ImGui::Button("Vector field setup")) {
	        vf->setup();
	        if (show_vector_field) display_vector_field();
	        update_mesh_colors();
	    }*/

	    ImGui::InputDouble("VF rotation angle (degrees)", &vf->angle_degrees);

	    if (ImGui::Button("Min energy (1 step)", ImVec2(205, 0))) {
	    	std::cout << "--- Minimizing energy, 1 step ---" << endl;
	    	if (!vf->is_setup) vf->setup();
	        VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	        cout << "Initial energy : " << vf->eval_E(X) << endl;
	        for (int ii = 0; ii < 1; ++ii){
	            vf->energy_single_step();
	            X = vf->squared ? vf->U.array().square() : vf->U;
				cout << "e = " << vf->eval_E(X) << endl;
	        }

	        X = vf->squared ? vf->U.array().square() : vf->U;
	        cout << "Final energy : " << vf->eval_E(X) << endl;
	        if (show_vector_field) display_vector_field();
	        update_mesh_colors();
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Button("(C) Curl reduction", ImVec2(-1, 0))) {
	        VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	        cout << endl << "Initial energy : " << vf->eval_E(X) << endl;
	        if (vf->is_setup) {
	            curl_reduction_one_step();
	        }
	        X = vf->squared ? vf->U.array().square() : vf->U;
	        cout << "Final energy : " << vf->eval_E(X) << endl << endl;
	    }

	    /*if (ImGui::Button("Min smoothness (UNconstrained) - 10 steps")) {
	        if (vf->is_setup) {
	            VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	            cout << "Initial smoothness : " << vf->eval_E_smoothness(X) << endl;
	            for (int ii = 0; ii < 10; ++ii)
	                vf->explicit_step_unconstrained(X);
	            cout << "Completed unconstrained smoothness step. Final smoothness : " << vf->eval_E_smoothness(X) << endl;

	            vf->get_U_from_X(X);
	            if (show_vector_field) display_vector_field();
	        }
	    }*/
	    if (ImGui::Button("(D) Vector field: Initial direct solve" , ImVec2(-1, 0))) {
	        if (!vf->is_setup) vf->setup();
	        vf->direct_solve();
	        if (show_vector_field) display_vector_field();
	        update_mesh_colors();
	    }

	    if (ImGui::Button("(V) Optimize vector field",ImVec2(-1, 0))) {
	        if (!vf->is_setup) vf->setup();
	        vf->minimize_energy();
	        if (show_vector_field) display_vector_field();
	        update_mesh_colors();
	    }
	    ImGui::Text("");

	    /*if (ImGui::Button("Normalize vector fields", ImVec2(205, 0)))
	        if (vf->is_setup)
	        {
	            vf->U.rowwise().normalize();
	            std::complex<double> complex_constant(3.0, 3.0);
	            vf->U = (vf->U.array() * complex_constant);
	            display_vector_field();
	        }
	    ImGui::SameLine(215);
	    if (ImGui::Button("Switch blue-green vector fields", ImVec2(-1, 0))) if (vf->is_setup) {
	        int nF = np3dp->mesh.F.rows(); assert(vf->U.size() == 2 * nF);
	        const VectorXcd v1 = vf->U.head(nF);
	        vf->U.head(nF) = vf->U.tail(nF);
	        vf->U.tail(nF) = v1;
	        display_vector_field();
	    }*/

	    ImGui::Text("--- Evaluate energies ---");
	    if (ImGui::Button("SmthVal", ImVec2(100, 0))) {
	        if (vf->is_setup) {
	            VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	            cout << "E_smoothness = " << vf->eval_E_smoothness(X) << endl;
	        }
	    }
	    ImGui::SameLine(110);
	    if (ImGui::Button("AlignVal", ImVec2(100, 0))) {
	        if (vf->is_setup) {
	            VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	            cout << "E_alignment = " << vf->eval_E_alignment(X) << endl;
	        }
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Button("OrthVal", ImVec2(100, 0))) {
	        if (vf->is_setup) {
	            VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	            cout << "E_orthogonality = " << vf->eval_E_orthogonality(X) << endl;
	        }
	    }
	    ImGui::SameLine(320);
	    /*if (ImGui::Button("CurlVal", ImVec2(100, 0))) {
	        if (vf->is_setup) {
	            double c = 0;
	            int nF = np3dp->F.rows();

	            VectorXi matching = vf->principal_matching(vf->U.head(nF), vf->effort);
	            c += vf->measure_curl_of_line_field(np3dp->EF, vf->U.head(nF), matching);
	            c += vf->measure_curl_of_line_field(np3dp->EF, vf->U.tail(nF), matching);
	            
	            cout << "Curl sum = " << c << endl;
	        }
	    }*/
	    if (ImGui::Button("UnitVal", ImVec2(100, 0))) {
	        if (vf->is_setup) {
	            VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;
	            cout << "E_unit = " << vf->eval_E_unit(X) << endl;
	        }
	    }

	    ImGui::Text("--- Check energy gradients using finite differences ---");
	    {
	        if (ImGui::Button("SmthGrad", ImVec2(100, 0))) {
	            if (vf->is_setup) {
	                VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;

	                VectorField* ptr_vf = vf.get();
	                std::function<double(const VectorXcd&)>    eval_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->eval_E_smoothness(X); };
	                std::function<VectorXcd(const VectorXcd&)> grad_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->dE_smoothness_dX(X); };

	                vf->check_analytic_gradient_using_finite_differences(X, eval_func, grad_func);
	            }
	        }
	        ImGui::SameLine(110);
	        if (ImGui::Button("AlignGrad", ImVec2(100, 0))) {
	            if (vf->is_setup) {
	                VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;

	                VectorField* ptr_vf = vf.get();
	                std::function<double(const VectorXcd&)>    eval_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->eval_E_alignment(X); };
	                std::function<VectorXcd(const VectorXcd&)> grad_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->dE_alignment_dX(X); };

	                vf->check_analytic_gradient_using_finite_differences(X, eval_func, grad_func);
	            }
	        }
	        ImGui::SameLine(215);
	        if (ImGui::Button("OrthGrad", ImVec2(100, 0))) {
	            if (vf->is_setup) {
	                VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;

	                VectorField* ptr_vf = vf.get();
	                std::function<double(const VectorXcd&)>    eval_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->eval_E_orthogonality(X); };
	                std::function<VectorXcd(const VectorXcd&)> grad_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->dE_orthogonality_dX(X); };

	                vf->check_analytic_gradient_using_finite_differences(X, eval_func, grad_func);
	            }
	        }
	        ImGui::SameLine(320);
	        if (ImGui::Button("UnitGrad", ImVec2(100, 0))) {
	            if (vf->is_setup) {
	                VectorXcd X = vf->squared ? vf->U.array().square() : vf->U;

	                VectorField* ptr_vf = vf.get();
	                std::function<double(const VectorXcd&)>    eval_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->eval_E_unit(X); };
	                std::function<VectorXcd(const VectorXcd&)> grad_func = [ptr_vf](const VectorXcd& X) { return ptr_vf->dE_unit_dX(X); };

	                vf->check_analytic_gradient_using_finite_differences(X, eval_func, grad_func);
	            }
	        }
	    }

	    ImGui::Text(" ");
	    if (ImGui::Button("Draw streamlines", ImVec2(200, 0))) {
	        prepare_viewer_index_data(vf_streamlines_viewer_index);

	        MatrixXd P1_u, P2_u, P1_v, P2_v;
	        if (vf->U.size() == 2 * np3dp->mesh.F.rows())
	            vf->draw_streamlines_of_U(P1_u, P2_u, P1_v, P2_v);

	        prepare_viewer_index_data(vf_streamlines_viewer_index);
	        viewer.data_list[vf_streamlines_viewer_index].line_width = 4.0;
	        viewer.data_list[vf_streamlines_viewer_index].add_edges(P1_u, P2_u, np3dp->dominant_strip_color * 0.9);
	        viewer.data_list[vf_streamlines_viewer_index].add_edges(P1_v, P2_v, np3dp->subdominant_strip_color * 0.9);
	    }
	    ImGui::SameLine(210);
	    if (ImGui::Button("Delete streamlines", ImVec2(200, 0))) {
			if (vf_streamlines_viewer_index >= 0) { viewer.data_list[vf_streamlines_viewer_index].clear();}
	    }

	    if (ImGui::BeginCombo("##combo", "Streamline tracing params")) // The second parameter is the label previewed before opening the combo.
		{
	        ImGui::InputInt("Streamline iterations", &vf->streamlines_iterations);
	        ImGui::InputDouble("Streamline dist ratio", &vf->streamlines_dist_ratio);
	        ImGui::InputDouble("Streamline dTime", &vf->streamlines_d_time);
		    ImGui::EndCombo();
		}

		if (ImGui::Button("Generate ROTATING streamlines", ImVec2(-1, 0))) {
	        prepare_viewer_index_data(vf_streamlines_viewer_index);
			cout << "Drawing rotating streamlines" << endl;
			vf->f_fraction = 0.9;
			vf->smoothness_coeff = 1;
			vf->alignment_coeff = 1;

			for (int i=0; i<181; ++i)
			{
				cout << " <<<<<<<<<<<<<<<<< a = " << i << endl;
				invalidate(true);
				vf->angle_degrees = i;

				vf->setup();
				vf->direct_solve();

		        MatrixXd P1_u, P2_u, P1_v, P2_v;
		        if (vf->U.size() == 2 * np3dp->mesh.F.rows())
		            vf->draw_streamlines_of_U(P1_u, P2_u, P1_v, P2_v);
			}
	    }


	    ImGui::Text(" ");
	    if (ImGui::Button("Optimize vf as cross-field", ImVec2(200, 0)))
	    {
	        if (!vf->is_setup) vf->setup();
	        vf->optimize_as_cross_field();
	        if (show_vector_field) display_vector_field();
	        update_mesh_colors();
	    }
		ImGui:ImGui::SameLine(210);
	    if (ImGui::Button("Integrate cross-field", ImVec2(200, 0))) {
	        if (vf->is_setup) {
	            np3dp->invalidate_parametrization_and_quad_related_data();

	            np3dp->UVcoords = vf->integrate_4_field_with_directional(true);
	            display_UV_coords(np3dp->meshCut.V, np3dp->meshCut.F, np3dp->UVcoords);

	            if (show_mesh_boundaries) display_mesh_boundaries(np3dp->meshCut.V, np3dp->meshCut.F);
	            np3dp->has_uv_coords = true;
	        }
	    }
	}

	////////////////////////////////////
	// --- Seamless Parametrization + Quad meshing
	{
	    ImGui::Text(" "); ImGui::Text("--- Seamless Parametrization + Quad meshing ---");

	    //ImGui::Checkbox("Fully seamless", &np3dp->fully_seamless_integration);
	    ImGui::InputDouble("Global scale", &np3dp->parametrization_global_scale);

	    if (ImGui::Button("Integrate 2x2-field", ImVec2(-1, 0))) {
	        if (vf->is_setup) {
	            np3dp->invalidate_parametrization_and_quad_related_data();

	            np3dp->UVcoords = vf->integrate_4_field_with_directional(false);
	            display_UV_coords(np3dp->meshCut.V, np3dp->meshCut.F, np3dp->UVcoords);

	            if (show_mesh_boundaries) display_mesh_boundaries(np3dp->meshCut.V, np3dp->meshCut.F);
	            np3dp->has_uv_coords = true;
	        }
	    }



	    if (ImGui::Button("Generate quad mesh", ImVec2(-1, 0))) {
	        if (np3dp->has_uv_coords) {
				np3dp->quad_mesh = QuadMesh(np3dp->output_folder);
	            np3dp->quad_mesh.generate_from_integration(np3dp->mesh, np3dp->meshCut, vf->intData, np3dp->UVcoords);
	            np3dp->has_quad_mesh = true;

	        	// find average length of edges in subdominant direction and set paths density based on it.
			    const double d = np3dp->quad_mesh.get_avg_len_edge_type(1) * np3dp->mesh_scale_factor;
			    np3dp->default_paths_density = static_cast<int>(round(d / np3dp->desired_avg_layer_height));
			    cout << "Average length of sub-dominant edges : " << d << ". Paths density : " << np3dp->default_paths_density << endl;

	            Helpers::write_Emap_to_txt(DATA_PATH + np3dp->output_folder + "Emap.txt", np3dp->quad_mesh.Emap, np3dp->quad_mesh.EV);

	            meshColors = FaceDegree;
	            show_quad_vertex_degrees = true;
	            set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);

	            update_mesh_colors(); display_quad_mesh_vertex_degrees(); display_edge_maps();
	            //display_UV_coords(np3dp->quad_mesh.V, np3dp->quad_mesh.T, np3dp->quad_mesh.UVcoords);

	            // remove unwanted visualizations
	            show_vector_field = false; if (vf_viewer_index >= 0)          viewer.data_list[vf_viewer_index].clear();
	            if (vf_streamlines_viewer_index >= 0)          viewer.data_list[vf_streamlines_viewer_index].clear();
	            if (vf_rotated_viewer_index >= 0)  viewer.data_list[vf_rotated_viewer_index].clear();
	            show_vf_sings = false; if (vf_sings_viewer_index >= 0)         viewer.data_list[vf_sings_viewer_index].clear();
	        }
	        else
	        {
	            cout << "The mesh does not have UV coords. First optimize the vector field, and then press on 'Integrate 2x2-field'. After that you can generate the quad mesh." << endl;
	        }
	    }
	}

	////////////////////////////////////
	// --- Ribs and directions
	if (np3dp->has_quad_mesh) 
	{
	    /* if (ImGui::Button("Go back to vector field stage", ImVec2(-1, 0))) {
	        np3dp->invalidate_parametrization_and_quad_related_data();
	        show_vector_field = true; show_vf_sings = true;
	        meshColors = MeshColors::VfConstraints;
	        if (show_vector_field) display_vector_field();
	        if (emap_viewer_index >= 0) viewer.data_list[emap_viewer_index].clear();
	        if (cuts_on_edges_viewer_index >= 0) viewer.data_list[cuts_on_edges_viewer_index].clear();
	        if (quad_vertex_degrees_viewer_index > 0)  viewer.data_list[quad_vertex_degrees_viewer_index].clear();
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        if (mark_fixed_pts_viewer_index >= 0) viewer.data_list[mark_fixed_pts_viewer_index].clear();
	        if (mark_helper_strips_viewer_index >= 0)  viewer.data_list[mark_helper_strips_viewer_index].clear();
	        if (alignment_of_selected_vertices_viewer_index >= 0)  viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
	        if (position_markers_viewer_index >= 0) viewer.data_list[position_markers_viewer_index].clear();
	        if (piece_labels_viewer_index >= 0) viewer.data_list[piece_labels_viewer_index].clear();
	        set_mesh_visual(np3dp->mesh.V, np3dp->mesh.F);
	        update_mesh_colors();
	    }
	    */

	    ImGui::Text(" "); ImGui::Text("--- Ribs and edge maps ---");

	    if (ImGui::Button("Find ribs positions", ImVec2(-1, 0))) {
	        np3dp->quad_mesh.cleanup();

	        // --- Create Eribs
			np3dp->quad_mesh.create_quad_ERibs();

	        // --- save results to txt
	        // Helpers::write_Emap_to_txt(DATA_PATH + np3dp->output_folder + "ERibs.txt", np3dp->quad_mesh.ERibMap, np3dp->quad_mesh.EV, false);
	        MatrixXd ribs = np3dp->quad_mesh.get_ribs_polylines();
			Helpers::write_matrix_to_txt(ribs, DATA_PATH + np3dp->output_folder + "/ribs_polylines.txt");
	        
	        show_ribs = true;
	        display_edge_maps();
	    }
	    if (ImGui::BeginCombo("##combo2", "Ribs params")) // The second parameter is the label previewed before opening the combo.
		{
	        ImGui::InputInt("Ribs_dist", &np3dp->quad_mesh.ribs_dist);
	        ImGui::InputInt("Ribs_f0_strip_network_U", &np3dp->quad_mesh.ribs_f0);
	        ImGui::InputInt("Ribs_f0_strip_network_V", &np3dp->quad_mesh.ribs_f0);
		    ImGui::EndCombo();
		}
	}


	////////////////////////////////////
	// --- Select side

	if (np3dp->has_quad_mesh)
	{
	    if (ImGui::Button("Flip quad Emap (green to blue and vice versa)", ImVec2(-1, 0)))
	    {
	        // flip emap
	        MeshHelpers::flip_Emap(np3dp->quad_mesh.Emap);

	        // flip UVcoords
	        VectorXd UV0 = np3dp->quad_mesh.UVcoords.col(0);
	        np3dp->quad_mesh.UVcoords.col(0) = np3dp->quad_mesh.UVcoords.col(1);
	        np3dp->quad_mesh.UVcoords.col(1) = UV0;

	        if (np3dp->has_partitioned_quad_mesh) display_partitioned_mesh();
	        display_edge_maps();

			np3dp->quad_mesh.strip_networks.invalidate();
	    	meshColors = FaceDegree;
			update_mesh_colors();
			if (strip_sequences_viewer_index >= 0) { viewer.data_list[strip_sequences_viewer_index].clear(); viewer.data_list[strip_sequences_viewer_index].clear_labels();}
	    }
	    ImGui::Text("");
	}

	ImGui::End();


	////////////////////////////////////
	// --- Mesh editing GUI
	ImGui::SetNextWindowPos(ImVec2(first_win_size+420, 0));
	ImGui::SetNextWindowSize(ImVec2(420, -1)); // 300
	ImGui::Begin("Quad Mesh editing");

	if (np3dp->has_quad_mesh)
	{
	    // ----------------------------------------------------------> Strips
	    ImGui::Text("Select:");
	    ImGui::SameLine(215);
	    ImGui::Checkbox("Forward rewire", &np3dp->quad_mesh.rewire_strip_forward);

		SelectionData& sel = np3dp->quad_mesh.selection;

	    if (ImGui::Checkbox("(S) Strip", &strip_selection_is_on)) {
	        if (strip_selection_is_on) {
	            vertex_selection_is_on = false; edge_selection_is_on = false; face_selection_is_on = false;
	            sel.selected_fi = -1;  sel.selected_ei = -1;
	        }
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
	        if (alignment_of_selected_vertices_viewer_index >= 0) viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
	        update_mesh_colors();
	    }
	    ImGui::SameLine(120);
	    if (ImGui::Button("Collapse (E)", ImVec2(90, 0))) {
	        if (sel.selected_strip.rows() == np3dp->quad_mesh.F.rows()) {
	            np3dp->quad_mesh.collapse_strip(sel.selected_strip, sel.selected_strip_direction);

	            np3dp->quad_mesh.strip_networks.invalidate();
	            update_visualization_after_topology_edits();
	        }
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Button("Subd (R)", ImVec2(90, 0))) {
	        if (sel.selected_strip.rows() == np3dp->quad_mesh.F.rows()) {
	            np3dp->quad_mesh.subdivide_strip(sel.selected_strip, sel.selected_strip_direction, true);

	            np3dp->quad_mesh.strip_networks.invalidate();
	            update_visualization_after_topology_edits();
	        }
	    }
	    ImGui::SameLine(310);
	    if (ImGui::Button("Rewire", ImVec2(90, 0))) {
	        if (sel.selected_strip.rows() == np3dp->quad_mesh.F.rows()){
	            np3dp->quad_mesh.rewire_strip(sel.selected_strip, sel.selected_strip_vis, sel.strip_is_closed, np3dp->quad_mesh.rewire_strip_forward);
	            np3dp->quad_mesh.strip_networks.invalidate();
	            update_visualization_after_topology_edits();
	        }
	    }

	    // ----------------------------------------------------------> Edges
	    if (ImGui::Checkbox("Edge", &edge_selection_is_on)) {
	        if (edge_selection_is_on) {
	            strip_selection_is_on = false; vertex_selection_is_on = false; face_selection_is_on = false;
	            sel.selected_fi = -1;  sel.selected_si = -1;
	        }
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        update_mesh_colors();
	    }
	    ImGui::SameLine(120);
	    if (ImGui::Button("Collapse", ImVec2(90, 0))) {
	        if (sel.selected_ei >= 0) {
	            np3dp->quad_mesh.collapse_edge(sel.selected_ei);

	            np3dp->quad_mesh.strip_networks.invalidate();
	            update_visualization_after_topology_edits();
	        }
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Button("Add Partition Edges", ImVec2(90, 0))) {
	        if (sel.selected_ei >= 0) {
	            vector<int>& eis = np3dp->manually_selected_partition_eis;
	            eis.push_back(sel.selected_ei);
	        }
	        draw_partitioning_cuts();
	    }
	    ImGui::SameLine(310);
		if (ImGui::Button("Clear partition edges", ImVec2(-1, 0))){
			 np3dp->manually_selected_partition_eis.clear(); draw_partitioning_cuts(); 
		}

	    // ImGui::Checkbox("Mark loop", &mark_edge_loop);
	    /*ImGui::SameLine(310);
	    if (ImGui::Button("Flatten edge loop", ImVec2(90, 0))) {
	        if (np3dp->selected_ei >= 0) 
	            QuadMeshEditing::flatten_edge_loop(np3dp.get(), np3dp->selected_ei);
	        np3dp->quad_mesh.has_strip_graphs = false;
	        update_visualization_after_topology_edits();
	    }*/

	    // ----------------------------------------------------------> Faces
	    if (ImGui::Checkbox("Face", &face_selection_is_on)) {
	        if (face_selection_is_on) {
	            strip_selection_is_on = false; vertex_selection_is_on = false; edge_selection_is_on = false; 
	            sel.selected_ei = -1; sel.selected_si = -1;
	        }
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        update_mesh_colors();
	    }
	    ImGui::SameLine(120);
	    if (ImGui::Button("Delete", ImVec2(90, 0))) {
			np3dp->quad_mesh.remember_previous_data();

	        // update Fquad and D
	        Helpers::removeRow(np3dp->quad_mesh.F, sel.selected_fi);
	        Helpers::removeRow(np3dp->quad_mesh.D, sel.selected_fi);

	        // cleanup, save and visualize
	        np3dp->quad_mesh.cleanup_unreferenced_vertices();
	        np3dp->quad_mesh.update_quad_mesh_data();
	        igl::writeOBJ(DATA_PATH + np3dp->output_folder + "quad_mesh.obj", np3dp->quad_mesh.V, np3dp->quad_mesh.F); cout << "Saved to file : " << DATA_PATH + np3dp->output_folder + "quad_mesh.obj" << endl;
	        np3dp->quad_mesh.create_quad_Emap();
	        Helpers::write_Emap_to_txt(DATA_PATH + np3dp->output_folder + "Emap.txt", np3dp->quad_mesh.Emap, np3dp->quad_mesh.EV);

	        np3dp->quad_mesh.strip_networks.invalidate();

	        update_visualization_after_topology_edits();
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Checkbox("4 vis", &sel.vis_4)){
	        sel.n_vis_to_select = sel.vis_4 ? 4 : 2;
	    }

	    ImGui::SameLine(310);
	    if (ImGui::Button("Add face", ImVec2(90, 0))) {
	        if (sel.selected_vis.size() > 3)
	        {
				np3dp->quad_mesh.remember_previous_data();

	            // update Fquad and D
	            int ii = np3dp->quad_mesh.F.rows();
	            np3dp->quad_mesh.F.conservativeResize(ii + 1, np3dp->quad_mesh.F.cols());
	            np3dp->quad_mesh.D.conservativeResize(ii + 1);
	            np3dp->quad_mesh.F.row(ii).setConstant(-1);
	            for (int i = 0; i < sel.selected_vis.size(); ++i)
	                np3dp->quad_mesh.F(ii, i) = sel.selected_vis[i];
	            np3dp->quad_mesh.D[ii] = 4;

	            // cleanup, save and visualize
	            np3dp->quad_mesh.cleanup_unreferenced_vertices();
	            np3dp->quad_mesh.update_quad_mesh_data();
	            igl::writeOBJ(DATA_PATH + np3dp->output_folder + "quad_mesh.obj", np3dp->quad_mesh.V, np3dp->quad_mesh.F); cout << "Saved to file : " << DATA_PATH + np3dp->output_folder + "quad_mesh.obj" << endl;
	            np3dp->quad_mesh.create_quad_Emap();
	            Helpers::write_Emap_to_txt(DATA_PATH + np3dp->output_folder + "Emap.txt", np3dp->quad_mesh.Emap, np3dp->quad_mesh.EV);

	            np3dp->quad_mesh.strip_networks.invalidate();

	            update_visualization_after_topology_edits();
	        }
	    }


	    // --- select vertices
	    if (ImGui::Checkbox("Vertices", &vertex_selection_is_on)) {
	        if (vertex_selection_is_on) {
	            strip_selection_is_on = false; edge_selection_is_on = false; face_selection_is_on = false;
	            sel.selected_fi = -1;  sel.selected_ei = -1; sel.selected_si = -1;
	        }
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
	    	if (mark_selected_vertices_viewer_index >= 0) { viewer.data_list[mark_selected_vertices_viewer_index].clear(); viewer.data_list[mark_selected_vertices_viewer_index].clear_labels(); }
			if (strip_sequences_viewer_index >= 0) { viewer.data_list[strip_sequences_viewer_index].clear(); viewer.data_list[strip_sequences_viewer_index].clear_labels(); }
	        update_mesh_colors();
	    }
	    ImGui::SameLine(120);
	    if (ImGui::Button("Link", ImVec2(90, 0))){
	        // If the vis are linked, when they become aligned, their alignment is 'fixed', i.e. no rewiring operations are permitted that would misalign them

	        if (sel.selected_vis.size() == 2){
	            cout << "Added vis : " << sel.selected_vis[0] << " , " << sel.selected_vis[1] << " for alignment tracking" << endl;
	            cout << "Once this pair becomes aligned, the alignment will be 'fixed', i.e. no rewiring operations are permitted that would misalign them" << endl;
	            np3dp->quad_mesh.vis_to_check_for_alignments.emplace_back(sel.selected_vis[0], sel.selected_vis[1]);

	            MatrixXd P1, P2; vector<int> vis;
	            for (const auto& a : np3dp->quad_mesh.aligned_vertices)
	                a.append_visualization_data(np3dp->quad_mesh.V, P1, P2, vis);
	            display_singularity_alignments(strip_sequences_viewer_index, P1, P2, vis);
	        }
	    }
	    ImGui::SameLine(215);
	    if (ImGui::Button("Clear links", ImVec2(90, 0))){
	        cout << "Cleared tracking of all alignments" << endl; np3dp->quad_mesh.vis_to_check_for_alignments.clear();
	        if (strip_sequences_viewer_index >= 0) viewer.data_list[strip_sequences_viewer_index].clear();
	    }

	    ImGui::SameLine(310);
	    if (ImGui::Button("Clear vis", ImVec2(-1, 0))) {
	        sel.selected_vis.clear();
	        update_visualization_after_topology_edits();
	    }

	    if (ImGui::Button("(Z) Undo", ImVec2(-1, 0))){
	        if (np3dp->quad_mesh.has_previous_data) {
				np3dp->quad_mesh.restore_previous_data();
	            np3dp->quad_mesh.strip_networks.invalidate();
	            update_visualization_after_topology_edits();
	        }
			else{
				cout << "No data to revert to" << endl;
			}
	    }
		// ImGui::SameLine(310); // subdivides only elongated strips
	    if (ImGui::Button("Subdivide elongated strips selectively", ImVec2(-1, 0))){
	        np3dp->quad_mesh.subdivide_whole_quad_mesh( true);
	        np3dp->quad_mesh.strip_networks.invalidate();
	        update_visualization_after_topology_edits();
	    }


	    //ImGui::SameLine(190);
	    if (ImGui::Button("Smoothen", ImVec2(100, 0))) {
	        np3dp->quad_mesh.smoothen(np3dp->mesh, np3dp->quad_mesh.smoothing_with_boundaries, np3dp->quad_mesh.smoothing_iterations, true);
			np3dp->quad_mesh.strip_networks.invalidate();
	        update_visualization_after_topology_edits();
	    }
	    ImGui::SameLine(110);
	    ImGui::Checkbox("With bndaries", &np3dp->quad_mesh.smoothing_with_boundaries);
		ImGui::SameLine(270);
		ImGui::Text("Iterations");
	    ImGui::SameLine(355);
	    ImGui::InputInt("Iterations", &np3dp->quad_mesh.smoothing_iterations);

	    if (ImGui::Button("Initialize strips graphs", ImVec2(290, 0))) {
	        Strips::create_strip_graphs(np3dp.get());

	        // --- singularity alignments
	        np3dp->quad_mesh.aligned_vertices.clear();
	        Strips::find_tracked_alignments(np3dp.get(), np3dp->quad_mesh.aligned_vertices);
			
	        // visualize
	        if (cuts_on_edges_viewer_index >= 0) viewer.data_list[cuts_on_edges_viewer_index].clear();
	        if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
	        if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
	     
	        MatrixXd P1, P2; vector<int> vis;
	        for (const auto& a : np3dp->quad_mesh.aligned_vertices)
	            a.append_visualization_data(np3dp->quad_mesh.V, P1, P2, vis);
	        display_singularity_alignments(strip_sequences_viewer_index, P1, P2, vis);
	     
	        // --- blocked strips
	        Strips::find_blocked_strips(np3dp.get(), np3dp->quad_mesh.aligned_vertices);
	        display_blocked_strips();

			if (save_graphviz)
			{
				Strips::write_graphviz(np3dp.get(), DATA_PATH + np3dp->output_folder + "Dgraph.dot", DATA_PATH + np3dp->output_folder + "_Dgraph.png", np3dp->quad_mesh.strip_networks.Dgraph);
				Strips::write_graphviz(np3dp.get(), DATA_PATH + np3dp->output_folder + "cDgraph.dot", DATA_PATH + np3dp->output_folder + "_cDgraph.png", np3dp->quad_mesh.strip_networks.cDgraph);
			}

	        vertex_selection_is_on = true;
	        strip_selection_is_on = false; edge_selection_is_on = false; face_selection_is_on = false;
	    }
		ImGui::SameLine(300);
		ImGui::Checkbox("Save graphviz", &save_graphviz);


	    if (np3dp->quad_mesh.strip_networks.has_strip_networks_data)
	    {
	        ImGui::Text(" "); ImGui::SameLine(100);

			SelectionData& sel = np3dp->quad_mesh.selection;

            if (ImGui::Button("Show connecting strips", ImVec2(200, 0))) {
                if (sel.selected_vis.size() == 2)
                {
                    int vi = sel.selected_vis[0];
                    int vj = sel.selected_vis[1];
                    vector<bool> aligned;
                    vector<vector<int>> alignment_mesh_paths;
                    Strips::find_strips_sequences_between_two_vertices(np3dp.get(), vi, vj, np3dp->quad_mesh.strip_networks.counterstrip_sequencies, np3dp->quad_mesh.strip_networks.vis_sequences, aligned, np3dp->quad_mesh.strip_networks.strip_sequencies, alignment_mesh_paths);

                    display_intermediate_strips_of_selected_vis(strip_sequences_viewer_index);
                }
            }
	        

	        ImGui::Text("Step >>>> "); ImGui::SameLine(100);
	        if (ImGui::Button("Decide", ImVec2(95, 0))) {
	            sel.selected_si = Strips::decide_on_next_step(np3dp.get());

	            if (sel.selected_si == -1) {
	                cout << "No next step was found" << endl;
	            }
	            else {
	                sel.selected_strip_direction = np3dp->quad_mesh.strip_networks.StoD[sel.selected_si];
	                const GraphVertexProps& strip_props = np3dp->quad_mesh.strip_networks.strip_index_to_properties(sel.selected_si);
	                sel.selected_strip = strip_props.strip;
	                sel.strip_is_closed = strip_props.strip_is_closed;
	                sel.selected_strip_vis = strip_props.strip_vis;
	           
	                // Display
	                if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
	                if (alignment_of_selected_vertices_viewer_index >= 0) viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
	                display_selected_strip_data();
	            }
	        }

	        ImGui::SameLine(200);
	        if (ImGui::Button("Do", ImVec2(100, 0))) {

				if (sel.selected_si != -1) // if a next step has been found
				{
		            if (mark_selected_viewer_index>=0) viewer.data_list[mark_selected_viewer_index].clear();

		            // perform selected operation
		            if (np3dp->quad_mesh.strip_networks.StoD[sel.selected_si] == 0)
		                np3dp->quad_mesh.collapse_strip(sel.selected_strip, sel.selected_strip_direction);
		            else
		                np3dp->quad_mesh.rewire_strip(sel.selected_strip, sel.selected_strip_vis, sel.strip_is_closed, np3dp->quad_mesh.rewire_strip_forward); // TODO!

		            // display
		            if (cuts_on_edges_viewer_index >= 0) draw_partitioning_cuts();
		            set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
		            display_edge_maps();
		            if (show_quad_vertex_degrees)  display_quad_mesh_vertex_degrees();

		            if (mark_selected_vertices_viewer_index >= 0){
		                viewer.data_list[mark_selected_vertices_viewer_index].clear();
		                for (int k=0; k<sel.selected_vis.size(); ++k)
		                {
							int vi = sel.selected_vis[k];
							viewer.data_list[mark_selected_vertices_viewer_index].add_points(np3dp->quad_mesh.V.row(vi), RowVector3d(1.0, 0.0, 0.0));
			                viewer.data_list[mark_selected_vertices_viewer_index].add_label(np3dp->quad_mesh.V.row(vi), std::to_string(k));
		                }
		                    
		            }

		            // --- invalidate strip graphs
		            if (mark_selected_viewer_index >= 0) viewer.data_list[mark_selected_viewer_index].clear();
		            if (mark_helper_strips_viewer_index >= 0) viewer.data_list[mark_helper_strips_viewer_index].clear();
		            if (strip_sequences_viewer_index >= 0) viewer.data_list[strip_sequences_viewer_index].clear();
		            if (alignment_of_selected_vertices_viewer_index >= 0) viewer.data_list[alignment_of_selected_vertices_viewer_index].clear();
	        		np3dp->quad_mesh.strip_networks.invalidate();
				}
				else
				{
					cout << "No step has been decided, no action will be done." << endl;
				}
	        }
	    }
	}
	ImGui::End();


	////////////////////////////////////
	// --- Partitioning (of both sides)

	ImGui::SetNextWindowPos(ImVec2(first_win_size + 2 * 420, 0));
	ImGui::SetNextWindowSize(ImVec2(420, -1)); // 300
	ImGui::Begin("Partitioning");

	if (np3dp->has_quad_mesh)
	{
	    {

	        //ImGui::InputInt("Number of connection rows", &num_of_connection_strip_rows);

	        // ImGui::SliderFloat("Proximity to wedge singularities", &np3dp->partition_proximity_to_wedge_singularities, 0.0f, 1.0f);

	        // ImGui::SliderAngle("+2*Pi Angle threshold snakes", &np3dp->snake_strips_angle_threshold);


	        if (ImGui::Button("Cuts", ImVec2(95, 0)))
	            draw_partitioning_cuts();
	        ImGui::SameLine(105);
	        if (ImGui::Button("Snakes", ImVec2(95, 0))){
	            VectorXd angle_sums;
	            Partitioning::identify_snake_strips(np3dp.get(), angle_sums);

	            for (int si = 0; si < angle_sums.size(); ++si)
	                if (angle_sums[si] < np3dp->snake_strips_angle_threshold)
	                    angle_sums[si] = np3dp->snake_strips_angle_threshold;

	            // face colors
			    int nF = np3dp->quad_mesh.T.rows();
			    assert(np3dp->quad_mesh.TF.rows() == nF);
			    MatrixXd face_colors = np3dp->mesh_base_color.replicate(nF, 1);

	            MatrixXd strip_colors;
				int count = 0;
				igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, angle_sums, true, strip_colors);
	            for (int si = 0; si < angle_sums.size(); ++si){
	            	if (angle_sums[si] <= np3dp->snake_strips_angle_threshold + 1e-16){
	            		strip_colors.row(si) = np3dp->mesh_base_color;
	            	}
					else{
						++count;
					}
	            }

			    for (int fi=0; fi<np3dp->quad_mesh.T.rows(); ++fi)
			    {
	                int f = np3dp->quad_mesh.TF[fi]; // quad face index

	                int si = np3dp->quad_mesh.strip_networks.FtoS(f, 0);
	                if (si >=0) {
	                    VectorXi strip = np3dp->quad_mesh.strip_networks.strips.row(si);
					    if (strip[f])
				            face_colors.row(fi) = strip_colors.row(si);
	                }
			    }
				cout << "Found " << count << " strips that might be snakes." << endl;
			    viewer.data().set_colors(face_colors);
	        }
	        ImGui::SameLine(210);
	        // if (ImGui::Button("Subdivide all strips", ImVec2(-1, 0))) {
	        //     np3dp->quad_mesh.subdivide_whole_quad_mesh(false);
	        //     np3dp->quad_mesh.strip_networks.invalidate();
	        //     np3dp->default_paths_density = ceil(np3dp->default_paths_density * 0.5); // reduce paths density to half
	        //     np3dp->quad_mesh.smoothen(np3dp->mesh, false, 1, true);
	        //     update_visualization_after_topology_edits();
	        //     np3dp->quad_mesh_is_subdivided += 1;
	        // }

	        if (ImGui::Button("Partition", ImVec2(-1, 0))) {
	            if (paths_viewer_index >= 0) viewer.data_list[paths_viewer_index].clear();
	            np3dp->quad_mesh.cleanup();

	            // first make sure that ERibMaps have been calculated
	            if (np3dp->quad_mesh.ERibMap.rows() != np3dp->quad_mesh.Emap.rows()) {// if Eribs have not been found, first find ribs
	                np3dp->quad_mesh.create_quad_ERibs();
	            }

	            // then partition each side + save resulting meshes
	            Partitioning::partition(np3dp.get());
	            np3dp->save_all_pieces_data(true);

	            // --- get alignment markings (if on the correct side)
	            //const int side_for_alignment_markings = np3dp->sides[0].nP() > np3dp->sides[1].nP() ? 0 : 1;
	            //Fabrication::set_alignment_markings(np3dp.get(), side_for_alignment_markings);

	            display_partitioned_mesh();
	        }

	    	ImGui::Checkbox("Stop at cut intersections", &np3dp->partitioning_stop_at_cut_intersections);
	        ImGui::SameLine(235);

	        
	        if (np3dp->has_partitioned_quad_mesh){
	        	if (ImGui::Button("Clear partitions", ImVec2(-1, 0)))
			        {
			            np3dp->pieces.clear();
			            np3dp->has_partitioned_quad_mesh = false;
						meshColors = FaceDegree;

        				set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
			            update_mesh_colors();
        				display_edge_maps();

			            if (cuts_on_edges_viewer_index >= 0) viewer.data_list[cuts_on_edges_viewer_index].clear();
			        }

		        if (ImGui::InputInt("Piece id", &np3dp->piece_ID)) {
					if (np3dp->piece_ID < 0 || np3dp->piece_ID >= np3dp->nP())
						np3dp->piece_ID = -1;

	                if (np3dp->piece_ID < 0 || np3dp->piece_ID >= np3dp->nP()) {
	                    set_mesh_visual(np3dp->quad_mesh.V, np3dp->quad_mesh.T);
	                    display_edge_maps();
						update_mesh_colors();
	                    // display_per_face_data_on_quad_mesh(np3dp->quad_mesh.T, np3dp->quad_mesh.TF, np3dp->Fmap);
	                    //display_UV_coords(np3dp->quad_mesh.V, np3dp->quad_mesh.T, np3dp->quad_mesh.UVcoords);
	                }
	                else {
	                    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
	                    hedra::polygonal_edge_topology(np3dp->pieces[np3dp->piece_ID].D, np3dp->pieces[np3dp->piece_ID].F, EV, FE, EF, EFi, FEs, innerEdges);
	                    MatrixXi T; VectorXi TF;
	                    hedra::triangulate_mesh(np3dp->pieces[np3dp->piece_ID].D, np3dp->pieces[np3dp->piece_ID].F, T, TF);
	                    set_mesh_visual(np3dp->pieces[np3dp->piece_ID].V, T);
	                    display_edge_maps();
	                    // display_per_face_data_on_quad_mesh(T, TF, np3dp->Ftypes[side][ID]);
	                    // display_UV_coords(np3dp->Vs[side][ID], np3dp->Ts[side][ID], np3dp->UVs[side][ID]);
	                }
	                display_paths();
		        }

	            if (ImGui::Button("Save strips of pieces", ImVec2(195, 0))){
	                for (int ID=0; ID < np3dp->pieces.size(); ++ID){
	                    const Piece& piece = np3dp->pieces[ID];

		                MatrixXi strips; // Strips x F (with values strips(si, fi) = 1/0 if fi belongs/doesn't belong to strip si)
						VectorXi strips_directions; // Strips x 1 (with values 0/1 for strip in the dominant/subdominant direction)
						vector<array<vector<int>, 2>> strips_vis;
						vector<bool> strips_are_closed;
						vector<vector<int>> eis, fis;

						Strips::get_all_strips(piece.V, piece.F, piece.D, piece.originalBndryFaces,
							piece.Emap, strips, strips_directions, strips_vis, strips_are_closed, eis, fis, MatrixXi(), vector<vector<set<int>>>());
						VectorXi strip_sums = strips.rowwise().sum();

						Helpers::write_matrix_to_txt(strips, Helpers::get_filename(DATA_PATH + np3dp->output_folder, "piece_strips", ID, "txt"));
						Helpers::write_vector_to_txt(strips_directions, Helpers::get_filename(DATA_PATH + np3dp->output_folder, "piece_strips_directions", ID, "txt"));
	                }
				}
	        }


	        // --- Path Tracing
	        if (np3dp->has_quad_mesh && np3dp->has_partitioned_quad_mesh)
	        {
	            {
	                ImGui::Text("--- Path tracing ---");
	                ImGui::InputInt("Paths density", &np3dp->default_paths_density);
	                if (ImGui::InputDouble("Desired avrg layer height", &np3dp->desired_avg_layer_height))
	                {
	                    // update default_paths_density
	                    double avg_length = 0;
	                    for (int ei = 0; ei < np3dp->quad_mesh.EV.rows(); ++ei) {
	                        avg_length += (np3dp->quad_mesh.V.row(np3dp->quad_mesh.EV(ei, 0)) - np3dp->quad_mesh.V.row(np3dp->quad_mesh.EV(ei, 1))).norm();
	                    }
	                    avg_length /= np3dp->quad_mesh.EV.rows();
	                    avg_length *= np3dp->mesh_scale_factor;
	                    np3dp->default_paths_density = avg_length / np3dp->desired_avg_layer_height;
	                }

	                // --- Path tracing on quad pieces
	                if (ImGui::Button("Trace paths", ImVec2(-1, 0))) {
	                    if (np3dp->has_partitioned_quad_mesh) {
	                        vector<int> Piece_IDs;
						    if (np3dp->piece_ID == -1 || np3dp->piece_ID >= np3dp->nP()){// then add all ids
						        for (int id = 0; id < np3dp->nP(); ++id)
						            Piece_IDs.push_back(id);
						    }
						    else{
							    Piece_IDs.push_back(np3dp->piece_ID);
						    }

	                        PathsTracing::trace_continuous_paths_on_quad_mesh(np3dp.get(), Piece_IDs);
	                        np3dp->has_paths_on_pieces = true;
	                        // display paths
	                        if (show_paths)
	                            display_paths();
	                    }
	                }
	            }
	        }


	    }
	}
	ImGui::End();

    ////////////////////////////////////

    ImGui::SetNextWindowPos(ImVec2(first_win_size + 3*420, 0));
    ImGui::SetNextWindowSize(ImVec2(420, -1)); // 300
    ImGui::Begin("Quad Mesh evaluation");
    if (np3dp->has_quad_mesh)
    {
        if (ImGui::Button("Edges length uniformity error", ImVec2(-1, 0))) {
            evaluate_edges_length_uniformity(0);
        }
        if (ImGui::Button("Orthogonality error", ImVec2(-1, 0))) {
            evaluate_orthogonality();
        }
        if (ImGui::Button("Alignment error", ImVec2(-1, 0))) {
            evaluate_alignment();
        }
        // ImGui::Text("");
        if (ImGui::Button("Quad planarity error", ImVec2(-1, 0)))
        {
            bool set_non_quads_to_zero = true;
            VectorXd deviation_from_planar = measure_quad_planarity_error(set_non_quads_to_zero);

            display_per_face_data_on_quad_mesh(np3dp->quad_mesh.T, np3dp->quad_mesh.TF, deviation_from_planar, true);
            Helpers::write_vector_to_txt(deviation_from_planar, DATA_PATH + np3dp->output_folder + "deviation_from_planar.txt");
        	cout << "Deviation from planarity. Max = " << deviation_from_planar.maxCoeff() << ", mean = " << deviation_from_planar.mean() << endl;

        }
    }


    ImGui::End();
    };
}



///////////////////////////////////////////////
/// --- Evaluate
///////////////////////////////////////////////

void Gui::evaluate_edges_length_uniformity(int dir)
{
    // --- get average edge length in each direction
    VectorXd lengths;

    for (int ei=0; ei<np3dp->quad_mesh.EV.rows(); ++ei){
	    if (np3dp->quad_mesh.Emap[ei] == dir){
            double d = (np3dp->quad_mesh.V.row(np3dp->quad_mesh.EV(ei, 0)) - np3dp->quad_mesh.V.row(np3dp->quad_mesh.EV(ei, 1))).norm();
            lengths.conservativeResize(lengths.size() + 1);
            lengths[lengths.size() - 1] = d;
	    }
    }

    double mean = lengths.mean();
    double variance = (lengths.array() - mean).square().sum() / (lengths.size() - 1);
    double standard_deviation = std::sqrt(variance) / mean;

    cout << "The standard deviations of the edge lengths are: " << standard_deviation << endl;
}

void Gui::evaluate_orthogonality()
{

    // --- sum up angle differences from Pi
    double sum = 0;
    int count = 0;

    for (int fi=0; fi<np3dp->quad_mesh.F.rows(); ++fi){ // --- in each face
        for (int k = 0; k < np3dp->quad_mesh.D[fi]; ++k){
	        int vi = np3dp->quad_mesh.F(fi, k);
            if (!np3dp->quad_mesh.boundaryVertices[vi]){ // --- for each non-boundary vertex
                int prev_vi = np3dp->quad_mesh.F(fi, (np3dp->quad_mesh.D[fi] + k - 1) % np3dp->quad_mesh.D[fi]);
                int next_vi = np3dp->quad_mesh.F(fi, (k + 1) % np3dp->quad_mesh.D[fi]);

                RowVector3d v1 = (np3dp->quad_mesh.V.row(prev_vi) - np3dp->quad_mesh.V.row(vi)).normalized();
                RowVector3d v2 = (np3dp->quad_mesh.V.row(next_vi) - np3dp->quad_mesh.V.row(vi)).normalized();

                double theta = acos(v1.dot(v2));
                sum += abs(0.5 * igl::PI - theta);
                ++count;
            }
        }
    }

    cout << "Deviation from orthogonality (0.5*PI) : average = " << sum / count << endl; // sum = " << sum << "
}

void Gui::evaluate_alignment()
{
    // This evaluation only considers constraints in the dominant (blue) direction
	if (!vf->is_setup)
	{
		cout << "The vector field is not setup! It must be setup to evaluate alignment. When the quad mesh is deserialized directly, the vector field is not setup to save time. Redo the vector field optimization and quad mesh generation from scratch to fix this." << endl;
		return;
	}

	// --- Project the vertices of quad mesh on original triangle faces
    VectorXd sqrD; // #P list of smallest squared distances
    VectorXi I; // #P list of primitive indices corresponding to smallest distances
    MatrixXd C; // #P by 3 list of closest points
    igl::point_mesh_squared_distance(np3dp->quad_mesh.V, np3dp->mesh.V, np3dp->mesh.F, sqrD, I, C);

    // --- For each edge of quad mesh, find the corresponding vector from original vector field and evaluate alignment 
    VectorXd alignments;
    VectorXd confidences;
    int count_dominant = 0;

    // --- create map of which faces are constrained
    VectorXi F_constrained; F_constrained.setZero(np3dp->mesh.F.rows());
    VectorXcd dirs; dirs.setZero(np3dp->mesh.F.rows());
    for (int i=0; i<vf->constraint_fis.rows(); ++i){
        int fi = vf->constraint_fis[i];
        if (fi < np3dp->mesh.F.rows()){ // only consider dominant direction
	        F_constrained[fi] = 1;
	        dirs[fi] = vf->constraint_dirs[i];
        }
    }

    for (int ei=0; ei<np3dp->quad_mesh.EV.rows(); ++ei)
    {
        if (np3dp->quad_mesh.Emap[ei] == 0) // only considering dominant direction
        {
            ++count_dominant;

        	int v0 = np3dp->quad_mesh.EV(ei, 0);
	        int v1 = np3dp->quad_mesh.EV(ei, 1);

	        RowVector3d proj_v0 = C.row(v0);
	        RowVector3d proj_v1 = C.row(v1);

	        int t0 = I[v0]; // triangle0
	        int t1 = I[v1]; // triangle1


        	if (F_constrained[t0] > 0 || F_constrained[t1] > 0) // only if the confidence weights are not 0
            {
        		RowVector3d b0 = np3dp->mesh.barycenters.row(t0); // barycenter of triangle 0
		        RowVector3d b1 = np3dp->mesh.barycenters.row(t1); // barycenter of triangle 1

		        double d0 = (proj_v0 - b0).norm(); // distance of projected v0 from barycenter0
		        double d1 = (proj_v1 - b1).norm(); // distance of projected v1 from barycenter1

        		double c0 = vf->confidence_weights[t0].real(); // confidence weights of two vectors
				double c1 = vf->confidence_weights[t1].real(); // confidence weights of two vectors

                Vector3d vec0, vec1; vec0.setZero(); vec1.setZero();

                if (F_constrained[t0] > 0)
					MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, vec0, dirs[t0], t0);
                if (F_constrained[t1] > 0)
        			MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, vec1, dirs[t1], t1);

                if (vec0.dot(vec1) < 0) // make sure vec0 and vec1 are pointing in the same direction
                    vec0 *= -1;

                Vector3d vf_vec = (d1 * vec0 + d0 * vec1) / (d0 + d1);
                double c = (d1 * c0 + d0 * c1) / (d0 + d1);

		        Vector3d edge_vec = static_cast<Vector3d>(proj_v0 - proj_v1);

        		alignments.conservativeResize(alignments.size() + 1);
                confidences.conservativeResize(confidences.size() + 1);

        		edge_vec.normalize();
                vf_vec.normalize();

		        alignments[alignments.size()-1] = abs(vf_vec.dot(edge_vec)); // Perfect alignment would be 0 (i.e. vector and edge completely orthogonal).
                confidences[confidences.size() - 1] = c;

                // cout << "edge_vec" << endl << edge_vec << endl;
                // cout << "vf_vec" << endl << vf_vec << endl;
                // cout << alignments[alignments.size() - 1] << endl;
                // cout << endl;
            }
        }
    }

    cout << "Deviation from alignment in constrained faces (in dominant dir only) : mean = " << (alignments*confidences).mean() << endl;
    cout << "(considered " << alignments.size() << " out of " << count_dominant << " dominant edges)" << endl;
}

VectorXd Gui::measure_quad_planarity_error(bool set_non_quads_to_zero)
{
    VectorXd planarity;
	hedra::planarity(np3dp->quad_mesh.V, np3dp->quad_mesh.D, np3dp->quad_mesh.F, planarity);

    if (set_non_quads_to_zero) {
        for (int fi = 0; fi < np3dp->quad_mesh.F.rows(); ++fi)
            if (np3dp->quad_mesh.D[fi] != 4)
                planarity[fi] = 0;
    }

    return planarity;
}
