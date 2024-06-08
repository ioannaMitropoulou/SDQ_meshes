//
//  Created by Ioanna Mitropoulou
//

#pragma once
#ifndef NP_3DP_PIECE_H
#define NP_3DP_PIECE_H

#include <Eigen/Core>

#include "paths.h"
#include "quad_mesh.h"
using namespace std;
using namespace Eigen;

// --- A Piece of the partitioned quad
class Piece
{
public:
    Piece(){} // constructor

	MatrixXd V; 
    MatrixXi F; 
    VectorXi D;
    MatrixXd faceNormals; // per face normals (consistent direction with the original mesh)
    VectorXi Emap; 
    VectorXi ERibMap; 
    VectorXi V_to_Vquad_original, F_to_Fquad_original; // keep track of which Vquad/Fquad each new vertex/face is coming from
	VectorXi originalBndryFaces; // Fquad x 1 : 0/1 if face isn't/is on the ORIGINAL boundary (i.e. of the uncut mesh)
	MatrixXd UV; // per vertex UVs
    MatrixXd GuidingDirs; // per face 'forward' direction

    VectorXi nP_per_strip; // Fx1 number of paths per strip (stored on every face of the strip)

	// --- Paths
    Paths::PathCollection pathsCollection;

    void get_avg_directions_and_extents(Vector3d& extents, Matrix3d& directions) const; // finds average forward (i.e. GuidingDir) and tangent direction, and the geometry extents along those directions
	MatrixXd get_per_face_normals(const MatrixXd& N_ref); // get per face normals (and make sure that they match the reference normals Nref of the original mesh) 
	MatrixXd get_guiding_directions();
	int piece_ei_to_original_mesh_ei(int piece_v0, int piece_v1, QuadMesh& quad_mesh) const; // return ei if edge was found, -1 otherwise

    void serialize(int piece_id, const std::string& serialize_folder) const;
	void deserialize(int piece_id, const std::string& serialize_folder);
};


#endif //NP_3DP_PIECE_H