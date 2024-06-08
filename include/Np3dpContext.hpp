//
//  Created by Ioanna Mitropoulou on 25.11.21.
//

#pragma once
#ifndef NP_3DP_CONTEXT_H
#define NP_3DP_CONTEXT_H

#include <Eigen/Core>
#include <directional/setup_integration.h>
#include <directional/readOBJ.h>

#include "piece.h"
#include "quad_mesh.h"


class Np3dpContext;
using namespace std;
using namespace Eigen;


// --- Np3dpContext: class that holds all data
class Np3dpContext {
public:
    Np3dpContext(const std::string& data_folder);

    bool found_correct_files = false;

    void invalidate_parametrization_and_quad_related_data();
    void serialize() const;
    void deserialize();


    // --------- Folders and data
    string data_folder;
    string output_folder;
    string serialize_folder;


    // --------- Triangle mesh
    directional::TriMesh mesh;

	double mesh_scale_factor = 0;
    double parametrization_global_scale = 0;


    // --------- Integration
    bool is_cut_and_combed = false;

    int texture_scale = 1;
    int texture_line_width = 6;

    bool has_uv_coords = false;
    directional::TriMesh meshCut; // mesh cut from singularities where seamless parametrization is applied
	MatrixXd UVcoords;


    // --------- Quad mesh data
    bool has_quad_mesh = false;
    QuadMesh quad_mesh;
    int quad_mesh_is_subdivided = 0;


    // --------- quad mesh partitioning
    bool has_partitioned_quad_mesh = false;

	VectorXi Fmap; // #Fquad x 1: stores the component ID (0, 1, ...) where each face belongs. -1 means that the face has not yet been assigned to a component

    vector<Piece> pieces;
    int piece_ID = -1;

	VectorXi PartitioningCuts; // Equad x 1 : 0/1 if the edge isn't/is a partitioning cut

	vector<int> manually_selected_partition_eis; // partition edges added by the user

    int nP() const { return pieces.size(); } // get number of pieces in side
    void save_all_pieces_data(bool quiet = false) const;

	float partition_proximity_to_wedge_singularities = 0.5;
	bool partitioning_stop_at_cut_intersections = false;

	float snake_strips_angle_threshold = igl::PI * 6;

	// --- Paths tracing
    int default_paths_density = 8; // density in each quad. Starting value that is updated
    bool has_paths_on_pieces = false;
    double desired_avg_layer_height = 2.0; // in mm


    // --- colors
	RowVector3d mesh_base_color = .95 * Eigen::RowVector3d::Ones(); // .89
	RowVector3d dominant_strip_color = RowVector3d(0.35, 0.50, 1.0);
	RowVector3d subdominant_strip_color = RowVector3d(0.35, 1.0, 0.50);
};

#endif //NP_3DP_CONTEXT_H