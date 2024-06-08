//
//  Created by Ioanna Mitropoulou on 11.12.21.
//

#pragma once
#include "Eigen/Core"
#include "Np3dpContext.hpp"
#include "paths.h"

using namespace std;
using namespace Eigen;


namespace PathsTracing
{
	void trace_continuous_paths_on_quad_mesh(Np3dpContext* np3dp, const vector<int>& Piece_IDs);

	void trace(Piece& p, int density);

	void assign_number_of_paths_per_strip(Piece& p, int default_paths_density, double mesh_scale_factor, double desired_avg_layer_height);


	// --- Single quad face
	vector<Paths::Segment> trace_segments_on_quad_face(int fi, int density, const MatrixXd& V, const MatrixXi& FH, const VectorXi& HE, const VectorXi& nextH,
														const VectorXi& HV, const VectorXi& Emap);

	vector<Paths::Segment> trace_segments_on_irregular_face(int fi, int density, const MatrixXd& V, const MatrixXi& F, const VectorXi& D, const VectorXi& Emap,
															const MatrixXd& UV, const MatrixXi& FH, const VectorXi& HE, const VectorXi& prevH, const MatrixXi& EV,
														    const RowVector3d& rot_guiding_dir, double avg_edge_len);

	// --- Helpers
	void cleanup(Paths::PathCollection& pathCollection);

	struct Intersections { // stores a pair of isoline zero-crossing on a face
		int e1 = -1, e2 = -1; // the indices of the intersected edges
		Vector3d p1, p2; // positions of intersections
		inline bool is_complete() const { return e1 != -1 && e2 != -1; }
		inline void add_intersection(int ei, const Vector3d& p) {
			if (e1 == -1) {
				e1 = ei; p1 = p;
			}
			else {
				e2 = ei; p2 = p;
			}
		}
	};

}
