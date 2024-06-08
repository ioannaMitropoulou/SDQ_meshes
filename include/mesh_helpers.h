//
// Created by Ioanna Mitropoulou on 29.11.21.
//

#pragma once
#ifndef NP_3DP_MESH_HELPERS_H
#define NP_3DP_MESH_HELPERS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <set>

#include "helpers.h"
#include <hedra/polygonal_write_OFF.h>
#include <igl/writeOBJ.h>

using namespace Eigen;
using namespace std;

namespace MeshHelpers {

	SparseMatrix<double> VtoF(const MatrixXd& V, const MatrixXi& F); // returns a FxV sparse matrix
	SparseMatrix<double> FtoV(const MatrixXd& V, const MatrixXi& F); // returns a VxF sparse matrix
	SparseMatrix<double> EtoF(const MatrixXi& F, const MatrixXi& EF); // returns a FxE sparse matrix

	int edge_index(int vi1, int vi2, const vector<vector<int>>& VE);
    bool halfedge_index(int v0, int v1, const vector<vector<int>>& VE, const MatrixXi& F, const MatrixXi& EF, int& fi, int& m);
    bool halfedge_index(int v0, int v1, int ei, const MatrixXi& F, const MatrixXi& EF, int& fi, int& m);

	void complex_local_to_cartesian_world_vf(const MatrixXd& B1, const MatrixXd& B2, MatrixXd& cartesian_vf, const VectorXcd& complex_vf);
    void complex_local_to_cartesian_world_vf(const MatrixXd& B1, const MatrixXd& B2, Vector3d& cartesian_vec, const complex<double>& complex_vec, int fi);
	void cartesian_world_to_complex_local_vf(const MatrixXd& B1, const MatrixXd& B2, const MatrixXd& cartesian_vf, VectorXcd& complex_vf);
    void cartesian_world_to_complex_local_vf(const MatrixXd& B1, const MatrixXd& B2, const Vector3d& cartesian_vec, complex<double>& complex_vec, int fi);

    Vector3d barycentric_coordinates(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& p);


	// --- cutting mesh
	void insert_path_to_halfedges_cut_matrix(const vector<int>& path, const MatrixXi& F, const MatrixXi& EF, const vector<vector<int>>& VE, const MatrixXi& TT,
        const MatrixXi& TTi, MatrixXi& C);
	void get_cuts_from_singularities(const MatrixXd& V, const MatrixXi& F, const MatrixXi& EF, const vector<int>& singularities_vf,
		const vector<vector<int>>& VE, set<int>& vertices_in_cut, MatrixXi& C);
	void get_cuts_edge_correspondence(const MatrixXd& V, const MatrixXi& F, const MatrixXi& EV, const MatrixXi& EF, const VectorXi& vertex_correspondence,
		const vector<vector<int>>& VE, VectorXi& ei_correspondence);


	// --- matching / checking halfedge quantities
    template <class SomeVector>
    inline void match_neighbor_he_quantities(const MatrixXi& EV, const MatrixXi& FE, const MatrixXi& EF, const std::string& description,
        bool alter_both_sides, SomeVector& he_vector) {
        // This function check the values for each edge's half-edges, and makes sure that they match.

        auto sign = [](double s)->int {return (0 < s) - (s < 0); };

        cout << " --- Matching halfedge quantities : " << description << endl;
        assert(he_vector.size() == 3 * FE.rows()); // the he_vector must be a quantity defined per he, i.e. 3 per face
        int count = 0;
        int count_sign_flips = 0;
        double total_magn = 0;

        for (int ei = 0; ei < EV.rows(); ++ei) {
            Vector2i edge_faces = EF.row(ei);
            int f1 = edge_faces[0];   int f2 = edge_faces[1];
            int epos_in_f1 = -1, epos_in_f2 = -1;
            if (f1 >= 0 && f2 >= 0) {

                if (FE(f1, 0) == ei) epos_in_f1 = 0;
                else if (FE(f1, 1) == ei) epos_in_f1 = 1;
                else if (FE(f1, 2) == ei) epos_in_f1 = 2;
                else assert(epos_in_f1 >= 0);

                if (FE(f2, 0) == ei) epos_in_f2 = 0;
                else if (FE(f2, 1) == ei) epos_in_f2 = 1;
                else if (FE(f2, 2) == ei) epos_in_f2 = 2;
                else assert(epos_in_f2 >= 0);

                double s1 = he_vector[3 * f1 + epos_in_f1];
                double s2 = he_vector[3 * f2 + epos_in_f2];

                int sign_s1 = sign(s1);
                int sign_s2 = sign(s2);

                // Now we make sure that these two are exactly the same value, even if in reality they are not due to
                // discrepancy of the projection from the left and from the right

                // --- magnitude matching
                if (abs(abs(s1) - abs(s2)) > 1e-8) { // if the two sides are not the same, get the avg
                    count += 1;
                    total_magn += abs(abs(s1) - abs(s2));

                    if (alter_both_sides) {
                        double avg_val = 0.5 * (abs(s1) + abs(s2));
                        he_vector[3 * f1 + epos_in_f1] = sign_s1 * avg_val; // sign * avg_val
                        he_vector[3 * f2 + epos_in_f2] = sign_s2 * avg_val; // sign * avg_val

                    }
                    else {
                        he_vector[3 * f1 + epos_in_f1] = sign_s1 * abs(s2); // sign * abs(value of f2)
                    }
                }

                // --- sign matching
                if (he_vector[3 * f1 + epos_in_f1] * he_vector[3 * f2 + epos_in_f2] > 0) {
                    // Then flip the value that is closer to zero
                    if (abs(s1) > (s2)) he_vector[3 * f2 + epos_in_f2] *= -1;
                    else                he_vector[3 * f1 + epos_in_f1] *= -1;
                    count_sign_flips += 1;
                }
                assert(he_vector[3 * f1 + epos_in_f1] * he_vector[3 * f2 + epos_in_f2] <= 0);

            }
        }
        cout << "Number of altered projections = " << count << " , Total magnitude = " << total_magn << endl;
        cout << "Number of sign flips = " << count_sign_flips << endl;
    }

    template <class SomeVector>
    inline void check_neighbour_he_quantities(const MatrixXi& EV, const MatrixXi& FE, const MatrixXi& EF, const std::string& description,
        SomeVector& he_vector) {

        cout << " --- Checking halfedge quantities : " << description << endl;
        assert(he_vector.size() == 3 * FE.rows()); // the he_vector must be a quantity defined per he, i.e. 3 per face

        for (int ei = 0; ei < EV.rows(); ++ei) {
            Vector2i edge_faces = EF.row(ei);
            int f1 = edge_faces[0];   int f2 = edge_faces[1];
            int epos_in_f1 = -1, epos_in_f2 = -1;
            if (f1 >= 0 && f2 >= 0) {
                if (FE(f1, 0) == ei) epos_in_f1 = 0;
                else if (FE(f1, 1) == ei) epos_in_f1 = 1;
                else if (FE(f1, 2) == ei) epos_in_f1 = 2;
                else assert(epos_in_f1 >= 0);

                if (FE(f2, 0) == ei) epos_in_f2 = 0;
                else if (FE(f2, 1) == ei) epos_in_f2 = 1;
                else if (FE(f2, 2) == ei) epos_in_f2 = 2;
                else assert(epos_in_f2 >= 0);

                double s1 = he_vector[3 * f1 + epos_in_f1];
                double s2 = he_vector[3 * f2 + epos_in_f2];

                if (abs(abs(s1) - abs(s2)) > 1e-6) {
                    cerr << " ----- Not matching sigmas! s1 = " << s1 << ", s2 = " << s2 << ". f1 = " << f1 << " , f2 = " << f2 << endl;
                    assert(false);
                }

            }
        }
        cout << "All helfedge quantities match." << endl;
    }

    // --- vertex groups collapsing
    void get_ground_boundary_islands_indices(const MatrixXd& V, const MatrixXi& F, const VectorXi& vs_to_original_boundary,
        const VectorXi& cutV_correspondence, VectorXi& ground_islands);

    void collapse_boundary(MatrixXd& Vc, MatrixXi& Fc, VectorXi vs_to_original_boundary_index, VectorXi cutV_correspondence,
        VectorXi& F_to_Fc, VectorXi& V_to_Vc);

    Vector3d face_gradient(const MatrixXd& V, const MatrixXi& F, int face_index, const Vector3d& face_normal,
        double face_area, const VectorXd& spinning_form); // not normalized gradient

    void flip_Emap(VectorXi& Emap);
    void flip_UVcoords(MatrixXd& UV);

    void get_boundary_polyline_and_vertices(MatrixXd& boundaryPolyline, VectorXi& boundary_vis, const MatrixXd& V, const MatrixXi& F);
    Vector3d projectToPolyline(const Vector3d& pt, const MatrixXd& Polyline);
    pair<Vector3d, bool> projectToBoundaryInDirection(const Vector3d& pt, const MatrixXd& boundaryPolyline, const Vector3d& dir);
}

#endif //NP_3DP_MESH_HELPERS_H