//
// Created by Ioanna Mitropoulou on 21.12.21.
//

#pragma once
#ifndef NP_3DP_PARTITIONING_H
#define NP_3DP_PARTITIONING_H

#include <Eigen/Core>
#include "Np3dpContext.hpp"

using namespace Eigen;
using namespace std;

namespace Partitioning {
	void partition(Np3dpContext* np3dp);

	void propagate_cuts_on_quad_mesh(Np3dpContext* np3dp, VectorXi& C);

	VectorXi create_Fmap(const VectorXi& D, const MatrixXi& F, const MatrixXi& EF, const MatrixXi& FE , const VectorXi& C);

	void separate_pieces(Np3dpContext* np3dp, const VectorXi& Emap_original, const MatrixXd& UV_original);

	void get_Emap_and_Erib_from_original_quad_mesh(VectorXi& Emap, VectorXi& Erib, const VectorXi& Emap_original, const VectorXi& ERibMap_original,
		const VectorXi& V_to_Vquad_original, const MatrixXi& EF, const MatrixXi& EV, const MatrixXd& Vquad_original, const MatrixXi& EVquad_original);
	// all 'original' data is coming from the unpartitioned mesh 

	void identify_snake_strips(Np3dpContext* np3dp, VectorXd& angle_sums);
}

#endif //NP_3DP_PARTITIONING_H