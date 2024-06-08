//
//  Created by Ioanna Mitropoulou on 25.11.21.
//

#pragma once
#ifndef NP_3DP_VECTOR_FIELD_H
#define NP_3DP_VECTOR_FIELD_H

#include "Eigen/Core"
#include <Eigen/Sparse>
#include "Np3dpContext.hpp"
#include <igl/opengl/glfw/Viewer.h>
#include "mesh_helpers.h"

using namespace std;
using namespace Eigen;

enum class VFGenerationType { CurvatureAligned = 0, BoundaryAligned, ConstraintsFromFile, Random };
enum class TimeSteps { Explicit = 0, Implicit };

class VectorField
{
public:
	VectorField(const shared_ptr<Np3dpContext>& np3dp): np3dp(np3dp), intData(2 * fieldDegree)
	{
		ftb.init(np3dp->mesh);
	}
	VectorField(const shared_ptr<Np3dpContext>& np3dp, TimeSteps time_steps, VFGenerationType vfGenerationType, double f_fraction, double smoothness_coeff, double alignment_coeff, double orthogonality_coeff, double unit_coeff, bool add_boundary_constraints, bool curvature_absolute_values): np3dp(np3dp),
		intData(2 * fieldDegree), time_steps(time_steps), vfGenerationType(vfGenerationType), f_fraction(f_fraction), smoothness_coeff(smoothness_coeff), alignment_coeff(alignment_coeff),
		orthogonality_coeff(orthogonality_coeff), unit_coeff(unit_coeff), add_boundary_constraints(add_boundary_constraints), curvature_absolute_values(curvature_absolute_values)
	{
		ftb.init(np3dp->mesh);
	}

	directional::IntrinsicFaceTangentBundle ftb;
	const bool squared = true;
	const int fieldDegree = 2;

	void setup();
	bool prepare_implicit_steps_solver();

    void save_vf_to_txt(const string& fileName, const VectorXcd& vf) const;


	// --- Parameters
	bool is_setup = false;
	int max_iterations = 50;
	TimeSteps time_steps = TimeSteps::Implicit ; //Implicit
	VFGenerationType vfGenerationType = VFGenerationType::CurvatureAligned;

	double decrease_mult = 0.8; // 1.0 means no decreasing.

	double smoothness_coeff = 10.0;
	double alignment_coeff = 0.1;
	double orthogonality_coeff = 2.0;
	double unit_coeff = 1.0;

	double f_fraction = 0.10;

	bool add_boundary_constraints = false;

	bool set_alignment_constraints_in_both_directions = true;
	bool set_smoothness_constraint_in_both_directions = true;
	bool set_unit_constraint_in_both_directions = true;

	// curvature-related params

	bool curvature_absolute_values = false;


	// --- Utils
	void serialize() const;
	void deserialize(bool do_setup);

	// --- Streamlines
	void VectorField::streamlines(const VectorXcd& vf, const std::string& output_filename, MatrixXd& P1, MatrixXd& P2);
	void draw_streamlines_of_U(MatrixXd& P1_u, MatrixXd& P2_u, MatrixXd& P1_v, MatrixXd& P2_v);
	int streamlines_iterations = 50;
	double streamlines_dist_ratio = 2.0;
	double streamlines_d_time = 0.5;

	// --- Solve : Update U and matching
	void minimize_energy();
	void direct_solve();
	void energy_single_step();

	// --- Curl
	double measure_U_curl(const VectorXi& matching);
	void project_U_curl_free(const VectorXi& matching);
	void principal_matching(const VectorXcd& vf); // fills in the matching and effort variables
	void curl_reduction_step(const VectorXi& matching, VectorXcd& vf) const;
	double measure_curl_of_line_field(const MatrixXi& EF, const VectorXcd& vf, const VectorXi& matching) const;

	// --- Visualization
	double sizeRatio = 0.2; //3.5;
	int sparsity = 0;
	double avgScale = 0.5;
	void display_vector_field(igl::opengl::glfw::Viewer& viewer, int vf_viewer_index, const VectorXcd& vf, int color_style) const;

	// --- Constraints
	VectorXi constraint_fis; // Nconstraints x 1 OR 2*Nconstraints x 1
	VectorXcd constraint_dirs; // Nconstraints x 1 OR 2*Nconstraints x 1
	VectorXcd confidence_weights; // Fx1 OR 2Fx1, real numbers (written as complex with complex part always 0). Remember: weights only apply for soft alignment

	void add_curvature_aligned_constraints();
	void add_constraints_on_naked_edges();
	void add_random_constraints();
	void add_constraints_from_file();

	// --- Utils
	double eval_E(const VectorXcd& X) const; // complete energy evaluation

	// --- E
	SparseMatrix<complex<double>> E_smoothness_matrix() const;
	SparseMatrix<complex<double>> E_alignment_matrix() const;
	SparseMatrix<complex<double>> E_orthogonality_matrix() const;
	SparseMatrix<complex<double>> E_unit_matrix() const;

	// --- dE/dX
	VectorXcd dE_smoothness_dX(const VectorXcd& X) const;
	VectorXcd dE_alignment_dX(const VectorXcd& X) const;
	VectorXcd dE_orthogonality_dX(const VectorXcd& X) const;
	VectorXcd dE_unit_dX(const VectorXcd& X) const;

	// --- eval E
	double eval_E_smoothness(const VectorXcd& X) const;
	double eval_E_alignment(const VectorXcd& X) const;
	double eval_E_orthogonality(const VectorXcd& X) const;
	double eval_E_unit(const VectorXcd& X) const;

	SparseMatrix<double> Constraints_matrix(int nF, const VectorXi& constraint_fis) const;


	// --- Time steps
	void explicit_step(VectorXcd& X) const;
	void implicit_step(VectorXcd& X) const;
	void explicit_step_unconstrained(VectorXcd& X) const; // step that ignores all alignment constraints

	// --- vector field to be found
	VectorXcd U; // Fx1 or 2Fx1 (i.e. np3dp.vf or [np3dp.vf and np3dp.vf_rotated])
	VectorXi matching; // Fx1
	VectorXd effort;

	void get_U_from_X(const VectorXcd& X); // updates U

	// --- Integration (pair of 2fields)
	directional::IntegrationData intData;
	directional::CartesianField get_2field_from_complex_vf(const VectorXcd& vf) const;
	void get_CCW_raw4field_from_U(MatrixXd& raw4Field) const;

	void setup_vf_integration(directional::TriMesh& meshCut, directional::CartesianField& combed_field, bool allow_all_cross_field_matchings);
	MatrixXd integrate_4_field_with_directional(bool allow_all_cross_field_matchings);

	void optimize_as_cross_field();

	void check_analytic_gradient_using_finite_differences(const VectorXcd& X, const std::function<double(const VectorXcd&)>& eval_func, const std::function<VectorXcd(const VectorXcd&)>& grad_func);

	// --- rotation
	double angle_degrees = 0;

private:
	shared_ptr<Np3dpContext> np3dp;

	double dt = 0.01; // Time steps

	const int curvature_confidence_type = 0; // 0: kmin-kmax, 1: sqrt(kmax^2 + kmin^2), 2: theta1 * (1 - exp(theta2 * (kmax-kmin)^2))
	const double s = 0.0; // 0.0;

	double n_s = 1.0; // normalize smoothness
	double n_a = 1.0; // normalize alignment
	double n_o = 1.0; // normalize orthogonality
	double n_u = 1.0; // normalize unit constraint
	double n_r = 1.0; // normalize regularizer

	bool prepare_curl_solver();
	void set_cols_of_matrix_to_zero(SparseMatrix<complex<double> >& M, int start_col, int end_col) const;

	SparseMatrix<double> Curl; // (E_int x 2F) real entries 
	
	SparseMatrix<complex<double>> A; // (E_int x F / 2E_int x 2F) : two entries on every row: 1, -r^2 (transfer to same coordinate system)
	SparseMatrix<double> W; // (E_int x E_int / 2E_int x 2E_int) - diagonal, per-edge weights: w = |e|^2 / 3*(A_f + A_g)
	SparseMatrix<complex<double>> K; // (FxF / 2Fx2F - diagonal) (has the weighted Gaussian curvature terms). Making it complex (with 0 as imaginary part), so that it can be multiplied with A without problems
	SparseMatrix<double> Mf, Mf_inv, Mf2_sqrt_inv; // (2Fx2F / 4Fx4F) - diagonal, per-face mass matrix
	SparseMatrix<complex<double>> Mfc; // 2Fx2F diagonal per-face mass matrix (but with complex entries)

	SparseMatrix<complex<double>> P1, P2; //P1: I | 0  , P2 : 0 | I

	// soft alignment
	VectorXcd R; // #Fx1 OR 2Fx1: soft alignment directions (squared)
	SparseMatrix<complex<double>> confidence_weights_mat;

	shared_ptr<SparseLU<SparseMatrix<complex<double>>>> implicit_steps_solver = nullptr; // SimplicialLDLT< SparseMatrix<complex<double> > >
	shared_ptr<LeastSquaresConjugateGradient<SparseMatrix<double>>> curl_solver = nullptr;
	
	const complex<double> J = complex<double>(0.0, 1.0);
};


#endif //NP_3DP_VECTOR_FIELD_H