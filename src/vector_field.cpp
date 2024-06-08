//
//  Created by Ioanna Mitropoulou on 25.11.21.
//

#include "vector_field.h"
#include <igl/gaussian_curvature.h>
#include <igl/eigs.h>
#include <igl/edge_topology.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>
#include <directional/directional_viewer.h>
#include <directional/principal_matching.h>
#include <directional/power_field.h>
#include <directional/streamlines.h>
#include <chrono>
#include <igl/principal_curvature.h>
#include <random>
#include <directional/integrate.h>
#include <directional/write_raw_field.h>
#include <directional/write_matching.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/false_barycentric_subdivision.h>
#include <directional/power_to_raw.h>

#include <directional/Deprecated/representative_to_raw.h>

void VectorField::setup()
{
	cout << "\nVector Field generation, setting up" << endl;

	MatrixXd& V = np3dp->mesh.V;
	MatrixXi& F = np3dp->mesh.F;
	MatrixXi& EF = np3dp->mesh.EF;
	MatrixXi& EV = np3dp->mesh.EV;
	MatrixXd& B1 = np3dp->mesh.FBx;
	MatrixXd& B2 = np3dp->mesh.FBy;

	int nF = F.rows();

	if (U.size() == 0)
	{
		U.setZero(2 * nF);
		U.head(nF) = VectorXcd::Constant(np3dp->mesh.F.rows(), 1, complex<double>(1, 0));
		U.tail(nF) = VectorXcd::Constant(np3dp->mesh.F.rows(), 1, complex<double>(0, 1));
	}

	// --- Initialize constraints if they are not already filled in
	{
		if (constraint_fis.size() == 0)
		{
			constraint_fis.setZero(0);
			constraint_dirs.setZero(0);

			confidence_weights.setZero(nF);

			switch (vfGenerationType) {

				case VFGenerationType::CurvatureAligned: {
					add_curvature_aligned_constraints(); break;
				}
				case VFGenerationType::BoundaryAligned: {
					add_constraints_on_naked_edges(); break;
				}
				case VFGenerationType::ConstraintsFromFile: {
					add_constraints_from_file(); break;
				}
				case VFGenerationType::Random: {
					add_random_constraints(); break;
				}
				default:
					break;
			}

			if (add_boundary_constraints) add_constraints_on_naked_edges();

			constexpr bool rotate_degrees = true;
			if (rotate_degrees)
			{
				double theta = angle_degrees * igl::PI / 180.0;
				complex<double> rotation_vec = complex<double>(cos(theta), sin(theta));

				for (int k=0; k<constraint_dirs.size(); ++k)
				{
					constraint_dirs[k] *= rotation_vec;
				}
			}

			// check and normalize constraints
			if (constraint_fis.rows() != constraint_dirs.rows()) { cerr << "Attention! At vf->setup, constraint_fis.rows() != constraint_dirs.rows()" << endl;}
			for (int i = 0; i < constraint_dirs.size(); ++i) {
				assert(norm(constraint_dirs[i]) > 0);
				constraint_dirs[i] /= sqrt(norm(constraint_dirs[i])); // normalize
			}

			// duplicate constraints for vf rotated also // That's not really needed since we have the orthogonality constraints
			int nC = constraint_fis.size();
			
			confidence_weights.conservativeResize(2 * nF);
			confidence_weights.tail(nF) = VectorXd::Zero(nF);

			if (set_alignment_constraints_in_both_directions){
				constraint_dirs.conservativeResize(2 * nC);
				constraint_fis.conservativeResize(2 * nC);

				constraint_fis.tail(nC) = VectorXi::Constant(nC, nF) + constraint_fis.head(nC);
				constraint_dirs.tail(nC) = complex<double>(0.0, 1.0) * constraint_dirs.head(nC); // all dirs rotated by 90
				confidence_weights.tail(nF) = confidence_weights.head(nF);
			}

			// set constraints on U
			for (int i = 0; i < constraint_fis.size(); ++i) {
				U[constraint_fis[i]] = constraint_dirs[i];
			}

			// set constrained entries of U
			for (int i = 0; i < nC; ++i) {
				int fi = constraint_fis[i];
				U[fi] = constraint_dirs[i];
			}
		}
	}

	// --- constraints: create variables mask (for hard alignment), or R and confidence weights (for soft alignment)
	{
		int nC = constraint_fis.size();
		
		R.setZero(U.size()); // alignment directions (squared) R 
		for (int i = 0; i < nC; ++i) {
			int fi = constraint_fis[i];
			R[fi] = squared ? pow(constraint_dirs[i], 2) : constraint_dirs[i];
			// confidence_weights[fi] = complex<double>(1.0, 0.0);
		}
		Helpers::diag(confidence_weights, confidence_weights_mat);

	}


	// --- Create coordinate-systems transfer matrix A (E_int x F OR 2*E_int x 2F, complex entries): with ei_conj()^2, -ej_conj()^2 on every row
	// --- Create diagonal matrix W (E_int x E_int OR 2E_int x 2E_int - real, diagonal) with per-edge weights w = |e|^2 / 3*(A_f + A_g)
	// --- Create curl matrix (E_int x 2F, real entries)
	{
		vector< Triplet<complex<double> > > ATriplets;
		vector< Triplet<double>> CurlRealTriplets;
		vector<double> e_weights;
		int count_e_int = 0;
		int total_e_int = 0;

		// find total count e_internal
		for (int ei = 0; ei < EF.rows(); ++ei) {
			int fi = EF(ei, 0);
			int fj = EF(ei, 1);

			if (fi == -1 || fj == -1) continue;  //boundary edge
			else ++total_e_int;
		}

		for (int ei = 0; ei < EF.rows(); ++ei) {
			int fi = EF(ei, 0);
			int fj = EF(ei, 1);

			if (fi == -1 || fj == -1) continue;  //boundary edge

			RowVector3d e = V.row(EV(ei, 1)) - V.row(EV(ei, 0));

			// Compute the complex representation of e on the local coordinate system of fi
			RowVector2d vef = Vector2d(e.dot(B1.row(fi)), e.dot(B2.row(fi))).normalized();
			complex<double> ef(vef(0), vef(1));

			// Compute the complex representation of e on the local coordinate system of fj
			Vector2d veg = Vector2d(e.dot(B1.row(fj)), e.dot(B2.row(fj))).normalized();
			complex<double> eg(veg(0), veg(1));
		
			if (squared)
			{
				// vf
				ATriplets.emplace_back(count_e_int, fi, pow(conj(ef), 2));
				ATriplets.emplace_back(count_e_int, fj, -pow(conj(eg), 2));
				// vf rotated
				ATriplets.emplace_back(total_e_int + count_e_int, nF + fi, pow(conj(ef), 2));
				ATriplets.emplace_back(total_e_int + count_e_int, nF + fj, -pow(conj(eg), 2));
			}
			else
			{
				// vf
				ATriplets.emplace_back(count_e_int, fi, conj(ef)); //-r^2
				ATriplets.emplace_back(count_e_int, fj, -conj(eg)); //-r^2
				// vf rotated
				ATriplets.emplace_back(total_e_int + count_e_int, nF + fi, conj(ef)); //-r^2
				ATriplets.emplace_back(total_e_int + count_e_int, nF + fj, -conj(eg)); //-r^2
			}


			// store curl matrix entries
			CurlRealTriplets.emplace_back(count_e_int, fi,   vef[0]);
			CurlRealTriplets.emplace_back(count_e_int, nF+fi,   vef[1]);

			CurlRealTriplets.emplace_back(count_e_int, fj,   -veg[0]);
			CurlRealTriplets.emplace_back(count_e_int, nF+fj,   -veg[1]);

			// store edge weight
			e_weights.push_back(3 * e.squaredNorm() / (np3dp->mesh.faceAreas(fi,0) + np3dp->mesh.faceAreas(fj,0)));

			++count_e_int; 
		}
		assert(count_e_int == total_e_int);


		// --- Create matrix A (2*E_int x 2F)
		A.conservativeResize(2 * total_e_int, 2 * nF);
		A.setFromTriplets(ATriplets.begin(), ATriplets.end());

		// --- Create curl matrix (E_int x 2F) -always these dimensions (it is applied twice when vfOptType == VFOptType::Double)
		Curl.conservativeResize(total_e_int, 2 * nF);
		Curl.setFromTriplets(CurlRealTriplets.begin(), CurlRealTriplets.end());

		bool success = prepare_curl_solver();
		if (!success) return;

		// --- Create per-edge weights matrix W (E_int x E_int - diagonal)
		VectorXd b = Eigen::Map<VectorXd>(e_weights.data(), e_weights.size()); // convert std vector to eigen vector
		b.conservativeResize(2*count_e_int);
		b.tail(count_e_int) = b.head(count_e_int);
		Helpers::diag(b, W);
	}

	// --- Create K matrix (FxF OR 2Fx2F - diagonal) (has the weighted Gaussian curvature terms)
	{
		SparseMatrix<double> VA; // voronoi area
		VectorXd GaussianCurv;
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, VA);
		VectorXd va = VA.diagonal();
		igl::gaussian_curvature(V, F, GaussianCurv);

		VectorXcd K_entries(nF);
		for (int fi=0; fi<nF; ++fi) {
			int vi0 = F(fi,0);
			int vi1 = F(fi, 1);
			int vi2 = F(fi, 2);
			K_entries[fi] = complex<double> (np3dp->mesh.faceAreas(fi,0) * (GaussianCurv[vi0] / va(vi0) + GaussianCurv[vi1] / va(vi1) + GaussianCurv[vi2] / va(vi2)) / 3.0, 0.0);
		}
		K_entries.conservativeResize(2*nF);
		K_entries.tail(nF) = K_entries.head(nF);
		Helpers::diag(K_entries, K);
	}

	// --- Create Mass matrices
	{
		VectorXd Ad(2 * nF);
		Ad.head(nF) = np3dp->mesh.faceAreas.col(0);
		Ad.tail(nF) = np3dp->mesh.faceAreas.col(0);
		Helpers::diag(Ad.cwiseSqrt().cwiseInverse(), Mf2_sqrt_inv);

		VectorXd A2(2*nF);
		A2.head(nF) = np3dp->mesh.faceAreas.col(0);
		A2.tail(nF) = np3dp->mesh.faceAreas.col(0);

		Helpers::diag(A2, Mf);
		Helpers::diag(A2.cwiseInverse(), Mf_inv);

		VectorXcd Ac2(2 * nF);
		for(int i=0; i<A2.size(); ++i) Ac2[i] = complex<double>(A2[i], 0);
		Helpers::diag(Ac2, Mfc);
	}

	// --- Calculate *dt* for explicit/implicit steps
	// here I'm only using the first FxF part of the matrix (even if it's actually 2Fx2F)
	{
		// Compute the smallest eigenvalue of the generalized eigen value problem: A u = s B u. In our case E*X = s * Mf * X
		MatrixXd EigVs; // eigenvectors
		VectorXd S; // eigenvalues
		MatrixXcd E_full = E_smoothness_matrix();
		MatrixXcd E = E_full.block(0,0, nF, nF); // only get the first FxF part of the matrix (even if we are optimizing the augmented field)

		// unroll E into a matrix with real entries (2Fx2F)
		SparseMatrix<double> E2f, M2f;
		vector < Triplet<double>>E2fTriplets;
		VectorXd A_2f(2*nF);
		for (int fi=0; fi<nF; ++fi)
		{
			E2fTriplets.emplace_back(2*fi,   2*fi,     E(fi, fi).real());
			E2fTriplets.emplace_back(2*fi,   2*fi+1, - E(fi, fi).imag());
			E2fTriplets.emplace_back(2*fi+1, 2*fi,     E(fi, fi).imag());
			E2fTriplets.emplace_back(2*fi+1, 2*fi+1,   E(fi, fi).real());

			A_2f[2*fi]   = np3dp->mesh.faceAreas(fi,0);
			A_2f[2*fi+1] = np3dp->mesh.faceAreas(fi,0);
		}
		Helpers::diag(A_2f, M2f);
		E2f.conservativeResize(2*nF, 2*nF);
		E2f.setFromTriplets(E2fTriplets.begin(), E2fTriplets.end());


		// --- solve using igl
		/*
		bool success = igl::eigs(E2f, M2f, 1, igl::EIGS_TYPE_SM, EigVs, S); // extract the first eigenvector and eigenvalue
		if (!success) {
			cerr << "Attention! igl::eigs did not converge. Setting default dt = " << dt << endl;
		}
		else {
			double smallest_eigen_val = S[0];
			dt = 1.0/smallest_eigen_val;
			cout << "Smallest eigenvalue = " << smallest_eigen_val << " , dt = " << dt << endl;
		}*/


		// --- solve using spectra
		/*// https://spectralib.org/doc/classspectra_1_1symgeigsshiftsolver_3_01optype_00_01boptype_00_01geigsmode_1_1shiftinvert_01_4
		// Construct matrix operation objects using the wrapper classes
		// We are going to solve the generalized eigenvalue problem
		//     A * x = lambda * B * x,
		// where A is symmetric and B is positive definite

		SparseMatrix<double>& A = E2f;
		SparseMatrix<double>& B = M2f;
		using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
		using BOpType = Spectra::SparseSymMatProd<double>;
		OpType op(A, B);
		BOpType Bop(B);

		//Helpers::write_to_matlab_format(A, DATA_PATH + np3dp->output_folder, "A");
		//Helpers::write_to_matlab_format(B, DATA_PATH + np3dp->output_folder, "B");

		// Construct generalized eigen solver object, seeking (three) generalized
		// eigenvalues that are closest to zero. This is equivalent to specifying
		// a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
		Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert> geigs(op, Bop, 1, 5, 1e-8); // op, Bop, 3, 6, -1e-6

		// Initialize and compute
		geigs.init();
		int nconv = geigs.compute(Spectra::SortRule::LargestMagn);

		// Retrieve results
		Eigen::VectorXd evalues;
		//Eigen::MatrixXd evecs;
		if (geigs.info() == Spectra::CompInfo::Successful)
		{
			evalues = geigs.eigenvalues();
			//evecs = geigs.eigenvectors();

			std::cout << "Number of converged generalized eigenvalues: " << nconv << std::endl;
			std::cout << "Generalized eigenvalues: " << evalues << std::endl;
			//std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;
		}
		else {
			cout << "Could not find smallest eigenvalue. Using default dt = " << dt << endl;
		}*/

		// --- solve using matlab
    	Engine* engine;

        MatrixXd EV; // Eigenvectors of the laplacian (w. mass matrix)
        MatrixXd ED; // Eigenvalues of the laplacian (w. mass matrix)

        // Launch MATLAB
        igl::matlab::mlinit(&engine);

        // Send Laplacian matrix to matlab
        igl::matlab::mlsetmatrix(&engine, "A", E2f);

        // Send mass matrix to matlab
        igl::matlab::mlsetmatrix(&engine, "B", M2f);

        // Extract the first 5 eigenvectors
        igl::matlab::mleval(&engine, "[EV,ED] = eigs(A,B,5,'smallestabs')");
        // Turn eigenvalue diagonal matrix into a vector
        igl::matlab::mleval(&engine, "ED=diag(ED)");

        // Retrieve the result
        igl::matlab::mlgetmatrix(&engine, "EV", EV);
        igl::matlab::mlgetmatrix(&engine, "ED", ED);
		// cout << "ED : " << ED << endl;

        for (int i=0; i<ED.size(); ++i)
        {
	        if (ED(i,0) > 1e-10)
	        {
                dt = 1 / ED(i, 0);
				// dt = 100.0 / ED(i, 0);
                cout << "Smallest eigenvalue = " << ED(i, 0) << " , step size dt = " << dt << endl;
                break;
	        }
        }
	}

	// --- Initialize Permutation matrices for joint opt of vf and vf_rotated
	{
		vector<Triplet<complex<int>>> P1_triplets, P2_triplets;
		P1.conservativeResize(nF, 2*nF);
		P2.conservativeResize(nF, 2*nF);
		for (int fi=0; fi<nF; ++fi)
		{
			P1_triplets.emplace_back(fi, fi, complex<double>(1, 0));
			P2_triplets.emplace_back(fi, nF+fi, complex<double>(1, 0));
		}
		P1.setFromTriplets(P1_triplets.begin(), P1_triplets.end());
		P2.setFromTriplets(P2_triplets.begin(), P2_triplets.end());
	}

	// --- get energy normalization values
	{
		VectorXcd ones; ones.setOnes(A.cols());

		n_s = 1 / (ones.adjoint() * E_smoothness_matrix() * ones).norm(); // normalize smoothness coefficient
		n_a = 1 / (ones.adjoint() * E_alignment_matrix() * ones).norm(); // normalize alignment coefficient
		n_o = 1 / (ones.adjoint() * E_orthogonality_matrix() * ones).norm(); // normalize orthogonality coefficient
		n_u = 1 / (ones.adjoint() * E_unit_matrix() * ones).norm(); // normalize unit constraint coefficient
		//cout << "normalize smoothness = " << n_s << " , normalize alignment = " << n_a << " , normalize orthogonality = " << n_o << endl;
	}

	// --- prepare implicit steps
	if (time_steps == TimeSteps::Implicit) {
		bool success = prepare_implicit_steps_solver();
		if (!success){
			cerr << "--- !!! --- Implicit solver decomposition failed --- !!! --- \nSwitching to *explicit* steps" << endl;
			time_steps = TimeSteps::Explicit;
		}
	}

	// --- setup: true
	is_setup = true;
}


void VectorField::save_vf_to_txt(const string& fileName, const VectorXcd& vf) const {
	MatrixXd vec_matrix(np3dp->mesh.F.rows(), 3);
    for (int fi=0; fi<np3dp->mesh.F.rows(); ++fi) 
        vec_matrix.row(fi) = vf[fi].real()*np3dp->mesh.FBx.row(fi) + vf[fi].imag() * np3dp->mesh.FBy.row(fi);
    Helpers::write_matrix_to_txt(vec_matrix, DATA_PATH + np3dp->output_folder + fileName);
}


bool VectorField::prepare_implicit_steps_solver()
{
	SparseMatrix<complex<double>> Id(U.size(), U.size());
	Id.setIdentity();
	SparseMatrix<complex<double>> A_implicit;
	A_implicit = Id + 2 * dt * Mf_inv * ( smoothness_coeff *    n_s * E_smoothness_matrix() + 
								          orthogonality_coeff * n_o * E_orthogonality_matrix() + 
										  unit_coeff *          n_u * E_unit_matrix() +
		                                  alignment_coeff     * n_a * Mf * confidence_weights_mat
		);


	// cout << "Solver computing A_implicit" << endl;
	implicit_steps_solver = make_shared<SparseLU<SparseMatrix<complex<double>>>>();
	implicit_steps_solver->compute(A_implicit);
	if (implicit_steps_solver->info() != Success) { // check if decomposition succeeded
		
	}
	return implicit_steps_solver->info() == Success;
}


void VectorField::get_U_from_X(const VectorXcd& X)
{
	if (!squared) { U = X; return; }

	if (U.size()%2 != 0) throw; // U should have an even number for size (2*nF)
	int nF = np3dp->mesh.F.rows();
	U = X.cwiseSqrt(); // initialize U as simple sqrt

	const MatrixXd& N = np3dp->mesh.faceNormals;

	// make sure rotated field fits matching (should be oriented in a CCW way)
	for (int fi = 0; fi < nF; ++fi) {
		complex<double> A_complex = U[fi]; // original field
		complex<double> B_complex = U[fi + nF]; // rotated field
	
		Vector3d A, b, d;
		MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, A, A_complex, fi);
		MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, b, B_complex, fi);
		d = b * -1;
	
		Vector3d n = N.row(fi);
	
		if (A.cross(b).dot(n) > A.cross(d).dot(n)) {
			// do nothing, the U[fi + nF] is already pointing in the correct direction
		}
		else {
			U[fi + nF] *= -1; // flip direction (to maintain CCW order)
		}
	}
}


bool VectorField::prepare_curl_solver() {
	curl_solver = make_shared<LeastSquaresConjugateGradient<SparseMatrix<double>>>();
	curl_solver->analyzePattern(Curl.transpose()* Curl);
	if (curl_solver->info() != Success) cerr << "--- !!! --- curl_solver analyzePattern() failed --- !!! ---" << endl;
	else							    cout << "curl_solver analyzePattern() succeeded" << endl;
	return curl_solver->info() == Success;
}


void VectorField::direct_solve()
{ // remember that here we don't consider the unit and regularizer energies (because these depend on the previous solution that doesn't exist yet)
	assert(is_setup);

	VectorXcd X = squared ? U.array().square() : U;

	cout << " --- Direct solving Energy" << endl;

	VectorXcd constraint_dirs_X = squared ? constraint_dirs.array().square() : constraint_dirs;

	bool success = true;
	SparseMatrix<complex<double>> B = Mf * confidence_weights_mat;
	SparseMatrix<complex<double>> lhs;
	lhs = 2 * (smoothness_coeff *    n_s * E_smoothness_matrix() +    // smoothness
		       alignment_coeff  *    n_a * B +                        // alignment
		       orthogonality_coeff * n_o * E_orthogonality_matrix()); // orthogonality

	SparseLU<SparseMatrix<complex<double>>> solver;
	solver.compute(lhs);
	if (solver.info() != Success) { cerr << "--- !!! --- Solver decomposition failed --- !!! --- " << endl; success = false;}

	VectorXcd rhs = alignment_coeff *  n_a * (B * R + B.adjoint() * R);

	X = solver.solve(rhs);
	if (solver.info() != Success) { cerr << "--- !!! --- Solving failed! --- !!! ---" << endl; success = false; }
	

	if (success) {
		cout << "Solver succeeded" << endl;
		X.rowwise().normalize();
		get_U_from_X(X); //U = squared ? X.cwiseSqrt() : X;

		principal_matching(U.head(np3dp->mesh.F.rows()));

	} else {
		cerr << " --- !!! --- Could not do direct solve of the energies  --- !!! --- " << endl;
		throw;
	}

	cout << "Energy value after direct solve = " << eval_E(X) << endl;
}


// --- E
SparseMatrix<complex<double>> VectorField::E_smoothness_matrix() const
{
	SparseMatrix<complex<double> > M = A.adjoint() * W * A - s * K;  // i.e. A in the notes

	if (!set_smoothness_constraint_in_both_directions){
		int nF = M.cols(); nF *= 0.5;
		set_cols_of_matrix_to_zero(M, nF, 2 * nF); // remove from subdominant
	}

	return M;
}
SparseMatrix<complex<double>> VectorField::E_alignment_matrix() const
{
	return Mf * confidence_weights_mat; // i.e. B in the notes
}
SparseMatrix<complex<double>> VectorField::E_orthogonality_matrix() const
{
	return (P1 + P2).adjoint() * (P1 + P2); // i.e. O in the notes
}
SparseMatrix<complex<double>> VectorField::E_unit_matrix() const
{
	SparseMatrix<complex<double>> M = Mfc;

	if (!set_unit_constraint_in_both_directions){
		int nF = M.cols(); nF *= 0.5;
		set_cols_of_matrix_to_zero(M, nF, 2 * nF); // remove from subdominant
	}

	return M;
}


// --- dE/dX
VectorXcd VectorField::dE_smoothness_dX(const VectorXcd& X) const
{
	return 2 * E_smoothness_matrix() * X;
}
VectorXcd VectorField::dE_orthogonality_dX(const VectorXcd& X) const
{
	return 2 * E_orthogonality_matrix() * X;
}
VectorXcd VectorField::dE_alignment_dX(const VectorXcd& X) const
{
	SparseMatrix<complex<double>> B = Mf * confidence_weights_mat;
	return 2 * B * X - B.adjoint() * R  - B.adjoint() * R;
}
VectorXcd VectorField::dE_unit_dX(const VectorXcd& X) const
{
	VectorXcd X_unit = X;
	X_unit.rowwise().normalize();
	SparseMatrix<complex<double>> Mfc = E_unit_matrix();
	// return 2 * Mfc * X - Mfc * X_unit - Mfc.adjoint() * X_unit; // This is identical with the expression below because Mfc = Mfc.adjoint() since Mfc is a diagonal matrix and all complex entries have 0 imaginary part 
	return 2 * Mfc * (X - X_unit);
}


void VectorField::check_analytic_gradient_using_finite_differences(const VectorXcd& X, const std::function<double(const VectorXcd&)>& eval_func, const std::function<VectorXcd(const VectorXcd&)>& grad_func)
{
	double tol = 1e-4;
	double eps = 1e-5;

	VectorXcd numeric_grad; numeric_grad.setZero(X.size());
	for (int i=0; i<X.size(); ++i)
	{
		VectorXcd dX; dX.setZero(X.size());

		// calculate real part of gradient at [i]
		dX[i] = complex<double>(eps, 0);
		double real = (eval_func(X+dX) - eval_func(X-dX)) / (2 * eps);

		// calculate imaginary part of gradient at [i]
		dX[i] = complex<double>(0, eps);
		double imag = (eval_func(X+dX) - eval_func(X-dX)) / (2 * eps);

		numeric_grad[i] = complex<double>(real, imag);
	}

	VectorXcd analytic_grad = grad_func(X);

	VectorXcd error = analytic_grad - numeric_grad;
	VectorXd error_norm = error.rowwise().norm();
	// cout << "error_norm : " << endl << error_norm << endl;
	cout << "error sums = " << error_norm.sum() << endl;
}


// --- eval E
double VectorField::eval_E(const VectorXcd& X) const
{
	double s = smoothness_coeff * n_s * eval_E_smoothness(X);
	double a = alignment_coeff * n_a * eval_E_alignment(X);
	double o = orthogonality_coeff * n_o * eval_E_orthogonality(X);
	double u = unit_coeff * n_u * eval_E_unit(X);
	return s + a + o + u;
}
double VectorField::eval_E_orthogonality(const VectorXcd& X) const
{
	return (X.adjoint() * E_orthogonality_matrix() * X).norm();
}
double VectorField::eval_E_smoothness(const VectorXcd& X) const
{
	return (X.adjoint() * E_smoothness_matrix() * X).norm();
}
double VectorField::eval_E_alignment(const VectorXcd& X) const
{
	return ((X - R).adjoint() * E_alignment_matrix() * (X - R)).norm();
}
double VectorField::eval_E_unit(const VectorXcd& X) const
{
	VectorXcd X_unit = X;
	X_unit.rowwise().normalize();
	return ((X - X_unit).adjoint() * E_unit_matrix() * (X - X_unit)).norm();
}


SparseMatrix<double> VectorField::Constraints_matrix(int nF, const VectorXi& constraint_fis) const
{
	int nC = constraint_fis.size(); // number of constraints

	vector<Triplet<double> > triplets;
	for (int ci=0; ci<nC; ++ci) {
		int fci = constraint_fis[ci];
		triplets.emplace_back(ci, fci, 1);
	}
	SparseMatrix<double> M_constraints(nC, nF);
	M_constraints.setFromTriplets(triplets.begin(), triplets.end());
	return M_constraints;
}


void VectorField::set_cols_of_matrix_to_zero(SparseMatrix<complex<double> >& M, int start_col, int end_col) const
{
	vector< Triplet<complex<double> > > triplets;
	
	for (int k=0; k<M.outerSize(); ++k){
		for (SparseMatrix<complex<double>>::InnerIterator it(M,k); it; ++it){
			auto& v = it.value();
			int row = it.row();   // row index
			int col = it.col();   // col index (here it is equal to k)
			int index = it.index(); // inner index, here it is equal to it.row()

			if (!(start_col < col && col < end_col))
				triplets.emplace_back(row, col, v);
		}
	}

	M.setFromTriplets(triplets.begin(), triplets.end());
}


void VectorField::explicit_step(VectorXcd& X) const
{
	assert(is_setup);
	VectorXcd g;

		g =  smoothness_coeff    * n_s * dE_smoothness_dX(X) + 
			 orthogonality_coeff * n_o * dE_orthogonality_dX(X) +
			 unit_coeff          * n_u * dE_unit_dX(X) +
		     alignment_coeff     * n_a * dE_alignment_dX(X);

	X = X - Mf_inv * g * dt;
}


void VectorField::explicit_step_unconstrained(VectorXcd& X) const
{
	assert(is_setup);
	VectorXcd g;
	g = smoothness_coeff * n_s * dE_smoothness_dX(X) + orthogonality_coeff * n_o * dE_orthogonality_dX(X);

	X = - Mf_inv * g * dt + X;
}


void VectorField::implicit_step(VectorXcd& X) const
{
	assert(is_setup);

	VectorXcd Xprev = X;
	VectorXcd Xprev_u = X.rowwise().normalized();

	// VectorXcd rhs = X - dt * Mf_inv * (unit_coeff * n_u * (E_unit_matrix() * Xprev_u + E_unit_matrix().adjoint() * Xprev_u) +
	// 								   regularizer_coeff * n_r * (Mfc * Xprev   + Mfc.adjoint() * Xprev));
	// remember Mfc = Mfc.adjoint(), so this expression is identical as the one below 

	SparseMatrix<complex<double>> B = Mf * confidence_weights_mat;
	VectorXcd rhs = X - dt * Mf_inv * 2 * (alignment_coeff * n_a * ( - B * R - B.adjoint() * R)
		                                   -unit_coeff * n_u * E_unit_matrix() * Xprev_u);

	X = implicit_steps_solver->solve(rhs);

	if (implicit_steps_solver->info() != Success) { // check if decomposition succeeded
		cerr << "--- !!! --- Solver decomposition failed! --- !!! ---" << endl;
		throw;
	}
}
 

void VectorField::minimize_energy()
{
	if (is_setup) {
		cout << "Minimizing vector field energy" << endl;
		cout << endl << " --- Iterative solving, number of iterations = " << max_iterations << endl;

		auto start = chrono::steady_clock::now();

		// --- find direct solution of smoothness energy (saved in U) - also updates matching
		// direct_solve(); 

		VectorXcd X = squared ? U.array().square() : U;
		double current_energy = eval_E(X);
		cout << "Starting energy value : " << current_energy << " , Starting curl sum : " << measure_U_curl(matching) << endl;

		double prev_value = 1e14; // very big number
		double e = 1e-4; // threshold number

		// save initial values of coefficients
		double smoothness_coeff_init = smoothness_coeff;
		double orthogonality_coeff_init = orthogonality_coeff;
		double unit_coeff_init = unit_coeff;

		for (int i = 0; i < max_iterations; ++i)
		{
			// ------------------------------------------------- Minimization step(s), using X

			// --- energy minimization step (updates X)
			if (time_steps == TimeSteps::Explicit) explicit_step(X);
			else			                       implicit_step(X);

			// --- update coefficients
			smoothness_coeff *= decrease_mult;
			unit_coeff *= decrease_mult;
			orthogonality_coeff *= decrease_mult;

			if (time_steps == TimeSteps::Implicit) prepare_implicit_steps_solver();
			// cout << "smoothness_coeff = " << smoothness_coeff << ", unit_coeff = " << unit_coeff << ", regularizer_coeff = " << regularizer_coeff << endl;


			// ------------------------------------------------- get U and update matching
			get_U_from_X(X);

			principal_matching(U.head(np3dp->mesh.F.rows())); // matching and effort calculated here, they should NOT be updated after the curl step

			// ------------------------------------------------- curl reduction, using U and matching
			project_U_curl_free(matching);

			X = squared ? U.array().square() : U;
			//if (i%5 == 0)
			{
				const double current_energy = eval_E(X);

				cout << "Iteration = " << i << " , Total energy = " << current_energy << endl; // No need to measure curl, it's always zero after the projection // " , Curl sum = " << measure_U_curl(matching) << endl;
				if (abs(current_energy - prev_value) > e) {
					prev_value = current_energy;
				}
				else
				{
					cout << "Energy has converged after " << i << " iteratons." << endl;
					auto end = chrono::steady_clock::now();
					cout << "Elapsed time in milliseconds: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
					// cout << "Elapsed time in seconds: "<< chrono::duration_cast<chrono::seconds>(end - start).count() << " sec";
					break;
				}

				if (i == max_iterations - 1)
				{
					cout << "Energy has NOT converged" << endl;
					auto end = chrono::steady_clock::now();
					cout << "Elapsed time in milliseconds: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
					break;
				}
					
			}
		}
		cout << "Final Energy = " << eval_E(X) << " , Final curl sum = " << measure_U_curl(matching) << endl; // No need to measure curl, it's always zero after the projection //


		// restore altered coefficients to original values
		smoothness_coeff = smoothness_coeff_init;
		unit_coeff = unit_coeff_init;
		orthogonality_coeff = orthogonality_coeff_init;
	}
}


void VectorField::optimize_as_cross_field()
{
	directional::CartesianField field, powerField, combed_field;
	int N = 4;
	int Nconstraints = static_cast<int>(constraint_fis.rows()) * 0.5;
	
	// --- fill in constFaces
	VectorXi constFaces = constraint_fis.block(0, 0, Nconstraints, 1);
	
	// --- fill in constVectors
	MatrixXd constVectors;
	constVectors.setZero(Nconstraints, 3);
	for (int i=0; i<Nconstraints; ++i)
	{
		Vector3d v;
		MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, v, constraint_dirs[i], constFaces[i]);
		constVectors.row(i) = static_cast<RowVector3d>(v);
	}
	//Helpers::write_vector_to_txt(constFaces, DATA_PATH + np3dp->output_folder + "constFaces.txt");
	//Helpers::write_matrix_to_txt(constVectors, DATA_PATH + np3dp->output_folder + "constVectors.txt");

	// --- optimize cross-field
	directional::power_field(ftb, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(), -1.0), N, powerField);
	directional::power_to_raw(powerField,N,field, true);
	directional::principal_matching(field);
	directional::combing(field, combed_field);

	// --- copy result to U
	int nF = np3dp->mesh.F.rows();
	U.setZero(nF * 2);
	for (int fi = 0; fi < nF; ++fi)
	{
		complex<double> c;
		Vector3d cartesian_vec = combed_field.extField.block(fi, 0, 1, 3).transpose();
		MeshHelpers::cartesian_world_to_complex_local_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesian_vec, c, fi);
		U[fi] = c;
		U[fi + nF] = c * J; // rotate by 90 degrees
	}
	U.rowwise().normalize();
}


void VectorField::energy_single_step()
{
	if (is_setup)
	{
		VectorXcd X = squared ? U.array().square() : U;
		
		//cout << endl << "Single step. Starting energy value = " << eval_E(X) << endl;
		if (time_steps == TimeSteps::Implicit) implicit_step(X);
		else					               explicit_step(X);

		//cout << "Final energy value = " << eval_E(X) << endl;

		get_U_from_X(X);  //U = squared ? X.cwiseSqrt() : X;

		principal_matching(U.head(np3dp->mesh.F.rows()));
	}
}


void VectorField::principal_matching(const VectorXcd& vf)
{
	MatrixXd cartesian_field, raw_field;
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesian_field, vf);
	directional::representative_to_raw(np3dp->mesh.faceNormals, cartesian_field, fieldDegree, raw_field);

	directional::CartesianField field(ftb);
	field.N = fieldDegree;
	field.fieldType = directional::fieldTypeEnum::RAW_FIELD;
	field.set_extrinsic_field(raw_field);
	directional::principal_matching(field);

	matching = field.matching;
	effort = field.effort;

	// Helpers::write_vector_to_txt(field.matching, DATA_PATH + np3dp->output_folder + "matching.txt");
}


double VectorField::measure_U_curl(const VectorXi& matching)
{
	if (matching.rows() == 0){cerr << "Matching has not been initialized" << endl; }

	const int nF = np3dp->mesh.F.rows();
	double c = 0;
	c += measure_curl_of_line_field(np3dp->mesh.EF, U.head(nF), matching);
	c += measure_curl_of_line_field(np3dp->mesh.EF, U.tail(nF), matching);
	return c;
}


void VectorField::project_U_curl_free(const VectorXi& matching)
{
	int nF = np3dp->mesh.F.rows();
	VectorXcd U1 = U.head(nF); curl_reduction_step(matching, U1); U.head(nF) = U1;
	VectorXcd U2 = U.tail(nF); curl_reduction_step(matching, U2); U.tail(nF) = U2;
}


void VectorField::curl_reduction_step(const VectorXi& matching, VectorXcd& vf) const // updates vf
{
	int nE = np3dp->mesh.EF.rows();
	int nF = np3dp->mesh.F.rows();

	if (matching.size() != nE) { cerr << "--- Attention! Given matching has the wrong size: " << matching.size() << " , with nE = " << nE  << endl; throw;}

	// --- initial u0
	VectorXd u0(2 * nF);
	for (int fi = 0; fi < nF; ++fi) {
		u0[fi] = vf[fi].real();
		u0[nF+fi] = vf[fi].imag();
	}

	// --- matching
	// flip entries of Curl on the one side when matching is not trivial
	SparseMatrix<double> Matching;
	vector<Triplet<double>> matching_triplets;
	int count_e_internal = 0;
	for (int ei = 0; ei < matching.size(); ++ei){
		int fi = np3dp->mesh.EF(ei, 0);
		int fj = np3dp->mesh.EF(ei, 1);
		if (fi == -1 || fj == -1) continue;  //boundary edge
		matching_triplets.emplace_back(count_e_internal, fi, 1.0);
		matching_triplets.emplace_back(count_e_internal, nF+fi, 1.0);
		matching_triplets.emplace_back(count_e_internal, fj, matching[ei] == 0 ? 1.0 : -1.0);
		matching_triplets.emplace_back(count_e_internal, nF+fj, matching[ei] == 0 ? 1.0 : -1.0);
		++count_e_internal;
	}
	Matching.conservativeResize(Curl.rows(), Curl.cols()); // (E_int x 2F)
	Matching.setFromTriplets(matching_triplets.begin(), matching_triplets.end());

	SparseMatrix<double> C = Curl.cwiseProduct(Matching);

	// check the rank of lhs_dense
	{
		//MatrixXd C_dense = MatrixXd(Curl);
		//FullPivLU<MatrixXcd> lu_decompC(C_dense);
		//auto rankC = lu_decompC.rank();
		//std::cout << "C dimensions  : " << C_dense.rows() << " , " << C_dense.cols() << endl;
		//std::cout << "C RANK  = " << rankC << endl;
	}


	// --- Solve system
	SparseMatrix<double> A = C * Mf2_sqrt_inv; // (E_int x 2F), matrix A in olga's pdf

	VectorXd b = -C * u0;
	VectorXd w;

	// --- solve with min-norm solution
	{
		//curl_solver.compute(A * A.transpose()); // lhs : (E_int x 2F) * (2F x E_int) -> (E_int x E_int)
		//if (curl_solver.info() != Success) {
		//	cerr << "--- !!! --- Curl_solver compute() failed --- !!! ---" << endl; return;
		//}
		//auto z = curl_solver.solve(b); // b -> rhs: (E_int x 2F) * (2F x 1) = (E_int x 1)
		//if (curl_solver.info() != Success) {
		//	cerr << "--- !!! --- Curl reduction solving failed --- !!! ---" << endl; return;
		//}
		//w = A.transpose() * z; // (2F x E_int) * (E_int x 1)
	}

	
	// --- solve with least squares conjugate gradient
	{
		curl_solver->factorize(A);
		w = curl_solver->solve(b);
		if (curl_solver->info() != Success) { cout << "solver for curl minimum norm solution failed!" << endl; return; }
	}

	const double error = (A * w - b).norm(); if (error > 1e-12) cerr << "Attention! error of curl minimum norm solution = " << error << endl;

	VectorXd v = Mf2_sqrt_inv * w; // (2F x 2F) * (2F x 1)

	VectorXd u = u0+v;
	for (int fi = 0; fi < vf.size(); ++fi){
		vf[fi] = complex<double>(u[fi], u[nF + fi]);

		// if vector has too large norm, then normalize
		double norm = sqrt(std::norm(vf[fi]));
		if (norm > 1e6)
			vf[fi] /= norm;
	}
}


double VectorField::measure_curl_of_line_field(const MatrixXi& EF, const VectorXcd& vf, const VectorXi& matching) const
{
	int nF = vf.size();

	// --- matching
	// flip entries of Curl on the one side when matching is not trivial
	SparseMatrix<double> Matching;
	vector<Triplet<double>> matching_triplets;
	int count_e_internal = 0;
	for (int ei = 0; ei < matching.size(); ++ei) {
		int fi = EF(ei, 0);
		int fj = EF(ei, 1);
		if (fi == -1 || fj == -1) continue;  //boundary edge
		matching_triplets.emplace_back(count_e_internal,   fi, 1.0);
		matching_triplets.emplace_back(count_e_internal, nF + fi, 1.0);
		matching_triplets.emplace_back(count_e_internal,   fj, matching[ei] == 0 ? 1.0 : -1.0);
		matching_triplets.emplace_back(count_e_internal, nF + fj, matching[ei] == 0 ? 1.0 : -1.0);
		++count_e_internal;
	}
	Matching.conservativeResize(Curl.rows(), Curl.cols()); // (E_int x 2F)
	Matching.setFromTriplets(matching_triplets.begin(), matching_triplets.end());
	SparseMatrix<double> C = Curl.cwiseProduct(Matching);

	VectorXd u_0(2 * nF);
	for (int fi = 0; fi < nF; ++fi) {
		u_0[fi]    = vf[fi].real();
		u_0[nF+fi] = vf[fi].imag();
	}

	const VectorXd per_edge_curl = C * u_0;

	return per_edge_curl.norm();
}


void VectorField::display_vector_field(igl::opengl::glfw::Viewer& viewer, int vf_viewer_index, const VectorXcd& vf, int color_style) const {
	cout << "Displaying per-face vector field" << endl;
	int nF = np3dp->mesh.F.rows();

	// convert vf to a cartesian field
	MatrixXd cartesian_field;
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesian_field, vf);

	MatrixXd display_field(nF, 3 * fieldDegree);
	directional::representative_to_raw(np3dp->mesh.V, np3dp->mesh.F, cartesian_field, fieldDegree, display_field); // fvf = representative. display_field = raw_field

	// get colors
	MatrixXd glyphColors (nF, 3 * fieldDegree);

	if (color_style == 0) {
		//RowVector3d col = RowVector3d(0, 0.9, 0.2);
		RowVector3d col = RowVector3d(0, 0.2, 0.9);
		glyphColors.block(0, 0, glyphColors.rows(), 3) = col.replicate(nF, 1);
		if (fieldDegree == 2)
			glyphColors.block(0, 3, glyphColors.rows(), 3) = col.replicate(nF, 1);
	}
	else {
		RowVector3d col = RowVector3d(0, 0.9, 0.2);
		//RowVector3d col = RowVector3d(0, 0.2, 0.9);
		glyphColors.block(0, 0, glyphColors.rows(), 3) = col.replicate(nF, 1);
		if (fieldDegree == 2)
			glyphColors.block(0, 3, glyphColors.rows(), 3) = col.replicate(nF, 1);
	}

	MatrixXd VField, CField;
	MatrixXi FField;

	directional::glyph_lines_mesh(np3dp->mesh.barycenters, np3dp->mesh.faceNormals, np3dp->mesh.EF, display_field, glyphColors, sizeRatio, avgScale, VField, FField, CField, sparsity, 0.2);
	viewer.data_list[vf_viewer_index].set_mesh(VField, FField);
	viewer.data_list[vf_viewer_index].set_colors(CField);
	viewer.data_list[vf_viewer_index].show_lines = false;
}


// --- Add constraints

void VectorField::add_curvature_aligned_constraints(){

	if (f_fraction > 0)
	{
		const MatrixXd& V = np3dp->mesh.V;
		const MatrixXi& F = np3dp->mesh.F;
		int numV = V.rows();
		int numF = F.rows();
		Eigen::MatrixXd VCBary; VCBary.setZero(numV+numF,3);
		Eigen::MatrixXi FCBary; FCBary.setZero(3*numF,3);

		igl::false_barycentric_subdivision(V, F, VCBary, FCBary); // use curvature values on the barycenter of each face

		VectorXd kmin, kmax;
		MatrixXd dmax,dmin;
		igl::principal_curvature(VCBary, FCBary, dmax, dmin, kmax, kmin, 5,true);

		MatrixXd dmax_ = dmax.bottomRows(numF);
		MatrixXd dmin_ = dmin.bottomRows(numF);
		dmax = dmax_;
		dmin = dmin_;
		VectorXd kmax_ = kmax.bottomRows(numF);
		VectorXd kmin_ = kmin.bottomRows(numF);
		kmax = kmax_;
		kmin = kmin_;

		dmax.rowwise().normalize();
		dmin.rowwise().normalize();

		if (curvature_absolute_values) // we only use dmax later on, so only update that (dmin is not used)
		{
			for (int fi = 0; fi < np3dp->mesh.F.rows(); ++fi) {
				if (abs(kmin[fi]) > abs(kmax[fi])){
					dmax.row(fi) = dmin.row(fi);
				}
			}
		}

		// project dmax on the plane of each face
		const MatrixXd& N = np3dp->mesh.faceNormals;
		for (int fi=0; fi<N.rows(); ++fi) {
			dmax.row(fi) -= N.row(fi) * N.row(fi).dot(dmax.row(fi)); // project on plane with normal N
		}
		dmax.rowwise().normalize();

		VectorXd confidence;
		if (curvature_confidence_type == 0) // --- confidence = kmax-kmin	
			confidence = (kmax - kmin).cwiseAbs();
		else if (curvature_confidence_type == 1) // --- confidence = sqrt(kmax^2 + kmin^2)
			confidence = (kmax.cwiseProduct(kmax) + kmin.cwiseProduct(kmin)).cwiseSqrt(); //sqrt(k1^2 + k2^2)
		else { // --- confidence = theta1 * (1 - exp(theta2 * (kmax-kmin)^2))
			double theta1 = 0.8;
		    double theta2 = -0.014;
		    confidence.setZero(kmax.rows());
		    for (int i=0; i< kmax.rows(); ++i){
	   			confidence[i] = theta1 * (1 - exp(theta2 * (kmax[i] - kmin[i]) * (kmax[i] - kmin[i])));
		    }
		}

		// ---> debug
		//Helpers::write_vector_to_txt(kmax, DATA_PATH + np3dp->output_folder + "kmax.txt");
		//Helpers::write_vector_to_txt(kmin, DATA_PATH + np3dp->output_folder + "kmin.txt");
		//Helpers::write_matrix_to_txt(dmax, DATA_PATH + np3dp->output_folder + "dmax.txt");
		//Helpers::write_matrix_to_txt(dmin, DATA_PATH + np3dp->output_folder + "dmin.txt");
		//Helpers::write_vector_to_txt(confidence, DATA_PATH + np3dp->output_folder + "confidence.txt");
		// <---

		auto is_boundary_face = [](int fi, const MatrixXi& FE, const MatrixXi& EF)->bool {
			for (int m = 0; m < 3; ++m) {
				int ei = FE(fi, m);
				if (EF(ei, 0) == -1 || EF(ei, 1) == -1)
					return true;
			}
			return false;
		};

		struct CurvatureOnFaces {
			CurvatureOnFaces(int fi, const Vector3d& vec, double magnitude) : fi(fi), vec(vec), confidence(magnitude) {}
			int fi;
			Vector3d vec;
			double confidence;
		};

		// sort based on confidence (from higher to lower)
		vector<CurvatureOnFaces> curvature_on_faces;
		for (int fi = 0; fi < np3dp->mesh.F.rows(); ++fi) {
			bool consider_face = !is_boundary_face(fi, np3dp->mesh.FE, np3dp->mesh.EF);

			if (consider_face) {
				Vector3d vec = dmax.row(fi);
				curvature_on_faces.emplace_back(fi, vec, confidence[fi]); // with 90 degrees rotation
			}
		}
		sort(curvature_on_faces.begin(), curvature_on_faces.end(), [](const CurvatureOnFaces& a, const CurvatureOnFaces& b) {return a.confidence > b.confidence; });

		auto get_w = [](double f_magn_diff) {
			//double k1 = 0.8;
			//double k2 = -0.014;
			//return k1 * (1 - exp(k2*pow(f_magn_diff,2)));

			// return Helpers::remap_unbound(f_magn_diff, min, max, 0.1, 1.0);
			return f_magn_diff;
		};

		// --- how many constraints to set
		int min_count = 5; // constrain at least 5 faces
		int n_constraints = std::max(min_count, static_cast<int>(round(np3dp->mesh.F.rows() * f_fraction))); // constrain at least f_fraction*100% of all faces
		cout << "Number of face constraints to be added = " << n_constraints << ", which is " << static_cast<double>(n_constraints)/ static_cast<double>(numF) * 100 << "% of all faces." <<  endl;

		double min_confidence = curvature_on_faces[std::min(n_constraints, static_cast<int>(curvature_on_faces.size())) - 1].confidence;
		double max_confidence = curvature_on_faces[0].confidence;

		cout << "Min confidence used = " << min_confidence << ", Max confidence used = " << max_confidence << endl;

		constraint_fis.setZero(n_constraints);
		constraint_dirs.setZero(n_constraints);

		// --- set constraints in complex representation
		for (int i = 0; i < std::min(n_constraints, static_cast<int>(curvature_on_faces.size())); ++i) {
			if (curvature_on_faces[i].confidence < 0.1) cerr << " --- Attention! Using a low-confidence curvature direction!" << endl;
			constraint_fis[i] = curvature_on_faces[i].fi;
			// cout << curvature_on_faces[i].fi << endl;
			RowVector3d vec = curvature_on_faces[i].vec;
			constraint_dirs[i] = complex<double>(vec.dot(np3dp->mesh.FBx.row(curvature_on_faces[i].fi)), vec.dot(np3dp->mesh.FBy.row(curvature_on_faces[i].fi))); // convert to complex representation
			confidence_weights[curvature_on_faces[i].fi] = complex<double>(get_w(curvature_on_faces[i].confidence), 0); // convert to complex representation
		}
	}
}


void VectorField::add_constraints_on_naked_edges() {
	// --- create boundary constraints
	const MatrixXd& V = np3dp->mesh.V;
	const MatrixXi& F = np3dp->mesh.F;
	const MatrixXd& N = np3dp->mesh.faceNormals;

	vector<vector<Index>> L; // vector of indices of boundary loops
	igl::boundary_loop(F, L); // L: list of loops where L[i] = ordered list of boundary g_vertices in loop i

	struct BoundaryEdge {
		int v1=-1, v2=-1, fid=-1;
		RowVector3d dir;
	};

	
	for (auto loop : L) {

		// --- find boundary edges
		vector<BoundaryEdge> b_edges;
		int no_constraints_per_boundary = static_cast<int>(loop.size());// - 1;
		for (int i = 0; i < no_constraints_per_boundary; ++i) {
			b_edges.push_back(BoundaryEdge());

			int v1 = loop[i];
			int v2 = loop[(i + 1) % no_constraints_per_boundary];

			b_edges.back().v1 = v1;
			b_edges.back().v2 = v2;

			int fid = -1;
			for (int i = 0; i < F.rows() && b_edges.back().fid < 0; ++i) {// iterate until fid where v1, v2 belong is found
				Vector3i f = F.row(i);
				if ((v1 == f[0] || v1 == f[1] || v1 == f[2]) && 
					(v2 == f[0] || v2 == f[1] || v2 == f[2]))
					fid = i; // found fid
			}
			assert(fid >= 0);
			b_edges.back().fid = fid;
			//cout << fid << endl;

			bool tangent_constraint = false;


			RowVector3d normal = N.row(fid);
			RowVector3d dir = (V.row(v1) - V.row(v2)).normalized();
			if (tangent_constraint)
				b_edges.back().dir = dir;
			else
				b_edges.back().dir = normal.cross(dir).normalized(); // with 90 degrees rotation

			if(b_edges.back().dir.norm() < 1e-6){cerr << "While creating boundary constraint on face : " << fid << " dir.norm() < 1e-6."<< endl; throw; }
		}


		// --- find constraint direction and add new constraint to b, bc
		for (int ii=0; ii<b_edges.size(); ++ii) {
			const auto& e = b_edges[ii];

			// make sure face hasn't been added to constraints yet
			bool face_already_passed = false;
			for (int j = 0; j < constraint_fis.size(); ++j) {
				if (e.fid == constraint_fis(j)) {
					face_already_passed = true; break;
				}
			}

			if (face_already_passed)
				continue; // do NOT overwrite direction in already-constrained faces

			int k = static_cast<int>(constraint_fis.rows());
			constraint_fis.conservativeResize(k + 1);
			constraint_dirs.conservativeResize(k + 1);


			constraint_fis(k) = e.fid; // append face index to b
			constraint_dirs[k] = complex<double>(e.dir.dot(np3dp->mesh.FBx.row(e.fid)), e.dir.dot(np3dp->mesh.FBy.row(e.fid)));
			confidence_weights[e.fid] = complex<double>(1.0, 0.0);
		}
	}
}


void VectorField::add_random_constraints() {
	int nF = np3dp->mesh.F.rows();
	int n_constraints = std::max(5, int(nF * f_fraction));
	cout << "Generating random vector field with n_constraints = " << n_constraints << endl;

	constraint_fis.setZero(n_constraints);  // face indices
	constraint_dirs.setZero(n_constraints);  // tangent vectors in world coordinates

	// --- generate random tangent vectors and save them as constraints
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_int_distribution<int> distr_int(0, nF - 1);
	std::uniform_real_distribution<double> distr_double(-1, 1);
	std::uniform_real_distribution<double> distr_angle(0, igl::PI * 2);

	vector<int> passed_fis;

	for (int i = 0; i < n_constraints; ++i) {
		int fi = distr_int(eng);
		while (Helpers::in_vector(passed_fis, fi)) { // get new fi until you find one that has not been used yet
			fi = distr_int(eng);
		}

		constraint_fis[i] = fi;
		passed_fis.push_back(fi);

		Vector3d N = np3dp->mesh.faceNormals.row(fi);
		Vector3d v1 = np3dp->mesh.V.row(np3dp->mesh.F(fi, 0)) - np3dp->mesh.V.row(np3dp->mesh.F(fi, 1)).normalized();
		Vector3d v2 = v1.cross(N);
		Vector3d tangent_vec = distr_double(eng) * v1 + distr_double(eng) * v2;
		tangent_vec.normalize();
		auto z = complex<double>(tangent_vec.dot(np3dp->mesh.FBx.row(fi)), tangent_vec.dot(np3dp->mesh.FBy.row(fi)));

		constraint_dirs[i] = z * complex<double>(cos(distr_angle(eng)), sin(distr_angle(eng)));
		confidence_weights[fi] = complex<double>(1.0, 0.0);
	}

	// --- save results to file
	/*ofstream data(fileName);
	for (int i = 0; i < n_constraints; ++i) {
		int fi = constraint_fis[i];
		Vector3d vec = constraint_dirs.row(i);
		data << fi << " " << vec[0] << " " << vec[1] << " " << vec[2] << endl;
	}
	data.close();*/
}


void VectorField::add_constraints_from_file() {
	MatrixXd dirs; VectorXd ws;
	Helpers::read_vf_constraints(DATA_PATH + np3dp->data_folder + "/fvf_constraints.txt", constraint_fis, dirs, ws);

	// convert dirs (real 3D coords) to constraint_dirs (complex vector in local coordinates)
	constraint_dirs.setZero(constraint_fis.size());
	for (int i = 0; i < constraint_fis.size(); ++i) {
		int fi = constraint_fis[i];
		//cout << fi << endl;
		constraint_dirs[i] = complex<double>(dirs.row(i).dot(np3dp->mesh.FBx.row(fi)), dirs.row(i).dot(np3dp->mesh.FBy.row(fi)));
		confidence_weights[fi] = complex<double>(ws[i], 0.0);
	}
}


directional::CartesianField VectorField::get_2field_from_complex_vf(const VectorXcd& vf) const
{
	int nF = np3dp->mesh.F.rows();

	// --- get raw_field for both vf and vf_rotated
	MatrixXd cartesianVF; // (Vx2) vector fields in global coordinates (the representative vector on each face)
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesianVF, vf);

	MatrixXd rawField; // (Vx2) vector fields
	directional::representative_to_raw(np3dp->mesh.faceNormals, cartesianVF, fieldDegree, rawField);

	directional::CartesianField field(ftb);
	field.N = 2;
	field.fieldType = directional::fieldTypeEnum::RAW_FIELD;
	field.set_extrinsic_field(rawField);

	return field;
}

void VectorField::get_CCW_raw4field_from_U(MatrixXd& raw4Field) const
{
	int nF = np3dp->mesh.F.rows();

	// --- get raw_field for both vf and vf_rotated
	MatrixXd cartesianVF, cartesianVFRotated; // (Vx3) vector fields in global coordinates (the representative vector on each face)
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesianVF, U.head(nF));
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesianVFRotated, U.tail(nF));

	cartesianVF *= np3dp->parametrization_global_scale;
	cartesianVFRotated *= np3dp->parametrization_global_scale;

	MatrixXd rawField, rawFieldRot; // (Vx2) vector fields
	directional::representative_to_raw(np3dp->mesh.faceNormals, cartesianVF, fieldDegree, rawField);
	directional::representative_to_raw(np3dp->mesh.faceNormals, cartesianVFRotated, fieldDegree, rawFieldRot);


	// --- create combined raw 4-field (in CCW order)
	raw4Field.setZero(nF, 3 * 4);

	// rawField
	raw4Field.block(0, 0, nF, 3) = rawField.block(0, 0, nF, 3);   // A
	raw4Field.block(0, 6, nF, 3) = rawField.block(0, 3, nF, 3);   // C

	// rawFieldRotated : select the correct order of the rotated field, so that the result is CCW
	for (int fi = 0; fi < nF; ++fi) {
		RowVector3d A = rawField.block(fi, 0, 1, 3);
		RowVector3d b = rawFieldRot.block(fi, 0, 1, 3);
		RowVector3d d = rawFieldRot.block(fi, 3, 1, 3);
		RowVector3d n = np3dp->mesh.faceNormals.row(fi);
		if (A.cross(b).dot(n) > A.cross(d).dot(n))
		{
			raw4Field.block(fi, 3, 1, 3) = b;
			raw4Field.block(fi, 9, 1, 3) = d;
		}
		else
		{
			raw4Field.block(fi, 3, 1, 3) = d;
			raw4Field.block(fi, 9, 1, 3) = b;
		}
	}
}


void VectorField::setup_vf_integration(directional::TriMesh& meshCut,  directional::CartesianField& combed_field, bool allow_all_cross_field_matchings)
{
	int nF = np3dp->mesh.F.rows();
	directional::CartesianField U_2field_raw = get_2field_from_complex_vf(U.head(nF));
	directional::CartesianField Urot_2field_raw = get_2field_from_complex_vf(U.tail(nF));
	directional::write_raw_field(DATA_PATH + np3dp->output_folder + "U_2field.rawfield", U_2field_raw);
	directional::write_raw_field(DATA_PATH + np3dp->output_folder + "Urot_2field.rawfield", Urot_2field_raw);

	intData = directional::IntegrationData(2 * fieldDegree); // clear integration data

	MatrixXd raw4Field;
	get_CCW_raw4field_from_U(raw4Field); // get combination of vf and vf_rotated as real raw fields in CCW order

	directional::CartesianField field(ftb);
	field.N = fieldDegree * 2;
	field.fieldType = directional::fieldTypeEnum::RAW_FIELD;
	field.set_extrinsic_field(raw4Field);

	directional::write_raw_field(DATA_PATH + np3dp->output_folder + "4_field.rawfield", field);
	Helpers::write_matrix_to_txt(raw4Field, DATA_PATH + np3dp->output_folder + "raw4Field.txt");
	Helpers::write_matrix_to_txt(raw4Field, DATA_PATH + np3dp->output_folder + "raw4Field.txt");
	
	// --- get matching, effort, and singularities
	if (!allow_all_cross_field_matchings) // only allows matchings: 0/2 -> pair of non-interchangeable 2fields
	{
		field.matching = matching * 2; // double the matching
		field.effort = effort; 
		directional::effort_to_indices(field);
	}
	else // typical principal matching, allows all matchings (i.e. 0/1/2/3) -> cross-field
	{
		directional::principal_matching(field);
	}
	Helpers::write_vector_to_txt(field.matching, DATA_PATH + np3dp->output_folder + "matching.txt");

	// --- Integrate
	std::cout << "Setting up integration for a " << fieldDegree * 2 << "-field." << std::endl;

	directional::setup_integration(field, intData, meshCut, combed_field);

	intData.verbose = true;
	intData.integralSeamless = true;
	intData.roundSeams = false;
}


MatrixXd VectorField::integrate_4_field_with_directional(bool allow_all_cross_field_matchings)
{
	cout << endl << " 4-field integration" << endl;
	
	directional::CartesianField combedField(ftb);
	setup_vf_integration(np3dp->meshCut, combedField, allow_all_cross_field_matchings);
	
	MatrixXd UV, cornerWholeUV;
	directional::integrate(combedField, intData, np3dp->meshCut, UV, cornerWholeUV);
	std::cout << "Completed integration!" << std::endl;
	
	//Extracting the UV from [U,V,-U, -V];
	MatrixXd UVout = UV.block(0, 0, UV.rows(), 2);
	
	Helpers::write_matrix_to_txt(UVout, DATA_PATH + np3dp->output_folder + "UVcoords.txt");
	igl::writeOBJ(DATA_PATH + np3dp->output_folder + "mesh_cut.obj", np3dp->meshCut.V, np3dp->meshCut.F);

	return UVout;
}


// --- Serialize, deserialize

void VectorField::serialize() const {
	igl::serialize(decrease_mult, DATA_PATH + np3dp->serialize_folder + "vf.decrease_mult.igl");

	igl::serialize(smoothness_coeff, DATA_PATH + np3dp->serialize_folder + "vf.smoothness_coeff.igl");
	igl::serialize(alignment_coeff, DATA_PATH + np3dp->serialize_folder + "vf.alignment_coeff.igl");
	igl::serialize(orthogonality_coeff, DATA_PATH + np3dp->serialize_folder + "vf.orthogonality_coeff.igl");
	igl::serialize(unit_coeff, DATA_PATH + np3dp->serialize_folder + "vf.unit_coeff.igl");
	igl::serialize(f_fraction, DATA_PATH + np3dp->serialize_folder + "vf.f_fraction.igl");
	igl::serialize(add_boundary_constraints, DATA_PATH + np3dp->serialize_folder + "vf.add_boundary_constraints.igl");
	igl::serialize(set_alignment_constraints_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_alignment_constraints_in_both_directions.igl");
	igl::serialize(set_smoothness_constraint_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_smoothness_constraint_in_both_directions.igl");
	igl::serialize(set_unit_constraint_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_unit_constraint_in_both_directions.igl");
	igl::serialize(curvature_absolute_values, DATA_PATH + np3dp->serialize_folder + "vf.curvature_absolute_values.igl");


	igl::serialize(max_iterations, DATA_PATH + np3dp->serialize_folder + "vf.max_iterations.igl");
	igl::serialize(time_steps, DATA_PATH + np3dp->serialize_folder + "vf.time_steps.igl");
	igl::serialize(vfGenerationType, DATA_PATH + np3dp->serialize_folder + "vf.vfGenerationType.igl");

	igl::serialize(constraint_fis, DATA_PATH + np3dp->serialize_folder + "vf.constraint_fis.igl");
	igl::serialize(constraint_dirs, DATA_PATH + np3dp->serialize_folder + "vf.constraint_dirs.igl");
	igl::serialize(confidence_weights, DATA_PATH + np3dp->serialize_folder + "vf.constraints_weights.igl");

	igl::serialize(U, DATA_PATH + np3dp->serialize_folder + "vf.U.igl");
	igl::serialize(matching, DATA_PATH + np3dp->serialize_folder + "vf.matching.igl");
	igl::serialize(effort, DATA_PATH + np3dp->serialize_folder + "vf.effort.igl");

	igl::serialize(angle_degrees, DATA_PATH + np3dp->serialize_folder + "vf.angle_degrees.igl");
}


void VectorField::deserialize(bool do_setup) {
	igl::deserialize(decrease_mult, DATA_PATH + np3dp->serialize_folder + "vf.decrease_mult.igl");

	igl::deserialize(smoothness_coeff, DATA_PATH + np3dp->serialize_folder + "vf.smoothness_coeff.igl");
	igl::deserialize(alignment_coeff, DATA_PATH + np3dp->serialize_folder + "vf.alignment_coeff.igl");
	igl::deserialize(orthogonality_coeff, DATA_PATH + np3dp->serialize_folder + "vf.orthogonality_coeff.igl");
	igl::deserialize(unit_coeff, DATA_PATH + np3dp->serialize_folder + "vf.unit_coeff.igl");
	igl::deserialize(f_fraction, DATA_PATH + np3dp->serialize_folder + "vf.f_fraction.igl");
	igl::deserialize(add_boundary_constraints, DATA_PATH + np3dp->serialize_folder + "vf.add_boundary_constraints.igl");
	igl::deserialize(set_alignment_constraints_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_alignment_constraints_in_both_directions.igl");
	igl::deserialize(set_smoothness_constraint_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_smoothness_constraint_in_both_directions.igl");
	igl::deserialize(set_unit_constraint_in_both_directions, DATA_PATH + np3dp->serialize_folder + "vf.set_unit_constraint_in_both_directions.igl");
	igl::deserialize(curvature_absolute_values, DATA_PATH + np3dp->serialize_folder + "vf.curvature_absolute_values.igl");


	igl::deserialize(max_iterations, DATA_PATH + np3dp->serialize_folder + "vf.max_iterations.igl");
	igl::deserialize(time_steps, DATA_PATH + np3dp->serialize_folder + "vf.time_steps.igl");
	igl::deserialize(vfGenerationType, DATA_PATH + np3dp->serialize_folder + "vf.vfGenerationType.igl");

	igl::deserialize(constraint_fis, DATA_PATH + np3dp->serialize_folder + "vf.constraint_fis.igl");
	igl::deserialize(constraint_dirs, DATA_PATH + np3dp->serialize_folder + "vf.constraint_dirs.igl");
	igl::deserialize(confidence_weights, DATA_PATH + np3dp->serialize_folder + "vf.constraints_weights.igl");

	igl::deserialize(U, DATA_PATH + np3dp->serialize_folder + "vf.U.igl");
	igl::deserialize(matching, DATA_PATH + np3dp->serialize_folder + "vf.matching.igl");
	igl::deserialize(effort, DATA_PATH + np3dp->serialize_folder + "vf.effort.igl");

	igl::deserialize(angle_degrees, DATA_PATH + np3dp->serialize_folder + "vf.angle_degrees.igl");

	is_setup = false;
	if (do_setup)
		setup();
}


void VectorField::streamlines(const VectorXcd& vf, const std::string& output_filename, MatrixXd& P1, MatrixXd& P2)
{
	int N = 2;
	MatrixXd cartesian_field, raw_field;
	MeshHelpers::complex_local_to_cartesian_world_vf(np3dp->mesh.FBx, np3dp->mesh.FBy, cartesian_field, vf);
	directional::representative_to_raw(np3dp->mesh.faceNormals, cartesian_field, N, raw_field);
	directional::CartesianField field(ftb);
	field.N = fieldDegree;
	field.fieldType = directional::fieldTypeEnum::RAW_FIELD;
	field.set_extrinsic_field(raw_field);
	directional::principal_matching(field);

	directional::StreamlineData data;
	directional::StreamlineState state;
	VectorXi seedFaces = VectorXi(); // if it's empty, then a distribution is set by directional
	seedFaces.setZero(1); seedFaces[0] = 2144; // corner doubly curved: 513; // catenoid_half: 917;
	directional::streamlines_init(field, seedFaces, streamlines_dist_ratio, data, state);

	for (int i=0; i<streamlines_iterations; ++i)
	{
		// cout << "Streamlines tracing, iteration = " << i << endl;
		directional::streamlines_next(data, state, streamlines_d_time);
	}

	const std::vector<Eigen::RowVector3d>& segStart = state.segStart;
	const std::vector<Eigen::RowVector3d>& segEnd = state.segEnd;
	const std::vector<int>& segOrigFace = state.segOrigFace;

	// cout << "segStart.size() : " << segStart.size() << endl;
	// cout << "segEnd.size() : " << segEnd.size() << endl;

	P1.setZero(segStart.size(), 3);
	for (int i=0; i<segStart.size(); ++i)
		P1.row(i) = segStart[i];

	P2.setZero(segEnd.size(), 3);
	for (int i=0; i<segEnd.size(); ++i)
		P2.row(i) = segEnd[i];

	// --- save results to txt file
	cout << "Saving to file : " << output_filename << " ... ";
	ofstream outputFile(output_filename);
	if(outputFile.is_open())
	{
		outputFile<< P1.rows() << endl;
		for (int i=0;i<P1.rows();i++)
			outputFile << P1(i, 0) << " " << P1(i, 1) << " " << P1(i, 2) << " " << 
						  P2(i, 0) << " " << P2(i, 1) << " " << P2(i, 2) <<  " " << segOrigFace[i] << endl;

		cout << "Completed" << endl;
	}
	else
	{
		cout << "Could not open " << output_filename << endl;
	}
}


void VectorField::draw_streamlines_of_U(MatrixXd& P1_u, MatrixXd& P2_u, MatrixXd& P1_v, MatrixXd& P2_v)
{
	// --- dominant (blue) direction u
	{
		VectorXcd vf = U.head(np3dp->mesh.F.rows());
		std::string output_filename = DATA_PATH + np3dp->output_folder + "streamlines_U_"+ std::to_string(int(angle_degrees)) + ".txt";
		// std::string output_filename = DATA_PATH + np3dp->output_folder + "streamlines_U.txt";
		streamlines(vf, output_filename, P1_u, P2_u);
	}

	// --- subdominant (green) direction v
	{
		VectorXcd vf = U.tail(np3dp->mesh.F.rows());
		std::string output_filename = DATA_PATH + np3dp->output_folder + "streamlines_V_"+ std::to_string(int(angle_degrees)) + ".txt";
		// std::string output_filename = DATA_PATH + np3dp->output_folder + "streamlines_U.txt";
		streamlines(vf, output_filename, P1_v, P2_v);
	}
}

