//
// Created by Ioanna Mitropoulou on 21.12.21.
//

#include "partitioning.h"

#include <directional/polygonal_edge_topology.h>
#include <hedra/dcel.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <igl/heat_geodesics.h>

#include "mesh_helpers.h"
#include "quad_mesh.h"
#include "strips.h"



void Partitioning::partition(Np3dpContext* np3dp)
{
	std::cout << "Partitioning" << endl;

	QuadMesh& m = np3dp->quad_mesh;

	VectorXi& PartitioningCuts = np3dp->PartitioningCuts;
	Partitioning::propagate_cuts_on_quad_mesh(np3dp, PartitioningCuts);

	Helpers::write_matrix_to_txt(m.get_lines_matrix_from_cuts(PartitioningCuts), DATA_PATH + np3dp->output_folder + "current_C.txt");

	// --- create Fmap and F_connections_map ONLY if it hasn't already been created
	np3dp->Fmap = create_Fmap(m.D, m.F, m.EF, m.FE, PartitioningCuts);
	Helpers::write_vector_to_txt(np3dp->Fmap, DATA_PATH + np3dp->output_folder + "Fmap.txt");

	// --- separate pieces
	Partitioning::separate_pieces(np3dp, m.Emap, m.UVcoords);

	// --- debug --->
	// np3dp->sides[current_side].save_all_pieces_data(np3dp, current_side, true);
	// <--- 

	// --- Calculate normals + guiding directions for each piece
	MatrixXd N_ref = QuadMesh::per_face_normals(np3dp->quad_mesh.V, np3dp->quad_mesh.T, np3dp->quad_mesh.F, np3dp->quad_mesh.TF); // per face normals of original quad
	for (int ID=0; ID < np3dp->pieces.size(); ++ID) {
		Piece& p = np3dp->pieces[ID];
		p.faceNormals = p.get_per_face_normals(N_ref);
		p.GuidingDirs = p.get_guiding_directions();
		Helpers::write_matrix_to_txt(p.GuidingDirs, Helpers::get_filename(DATA_PATH + np3dp->output_folder, "guiding_dirs", ID, "txt"), true);
	}

	np3dp->has_partitioned_quad_mesh = true;
	std::cout << endl << "Completed partitioning with number of partitioned pieces = " << np3dp->nP() << endl << endl;
}


void Partitioning::propagate_cuts_on_quad_mesh(Np3dpContext* np3dp, VectorXi& C)
{
	bool stop_at_cut_intersections = np3dp->partitioning_stop_at_cut_intersections;
	float proximity_to_wedge_singularities = np3dp->partition_proximity_to_wedge_singularities;

	const QuadMesh& m = np3dp->quad_mesh;

	C.setZero(m.EV.rows()); // cuts (#E x 1) 0: not a cut, 1: cut

	vector<int> passed_vis;

	///////////////////////////////////////////
	// --- create cuts from irregular vertices

	deque<int> start_halfedges;
	for (int vi = 0; vi<m.V.rows(); ++vi)
	{
		if (m.vertexDegree[vi] != 4 && m.boundaryVertices[vi] != 1) { // if it's not a regular vertex, and it's not on the boundary
			cout << "Splitting from vi = " << vi << " , with vertex degree : " << m.vertexDegree[vi] << endl;
			
			if (m.vertexDegree[vi] == 2) { // singularity index = 1

				/*int he1, he2;
				if (Emap[HE[VH[vi]]] == 0)
				{
					he1 = VH[vi];
					he2 = twinH[nextH[nextH[nextH[he1]]]];
				}
				else
				{
					he2 = VH[vi];
					he1 = twinH[nextH[nextH[nextH[he2]]]];
				}
				start_halfedges.push_front(he1);
				// --- find other two starting cuts
				// first find how long
				int remember_he2 = he2;
				int d = QuadMesh::find_distance_to_boundary(he2, HV, HE, nextH, twinH);
				for (int i=0; i<d*proximity_to_wedge_singularities; ++i)
				{
					he2 = twinH[nextH[he2]]; if (he2 == -1) break;
					he2 = nextH[he2];
				}
				if (he2 >=0)
				{
					if (nextH[he2] >= 0)
						start_halfedges.push_front(nextH[he2]);
					if (twinH[prevH[twinH[he2]]] >= 0)
						start_halfedges.push_front(twinH[prevH[twinH[he2]]]);
				}
				start_halfedges.push_back(remember_he2);*/

				int he = m.VH[vi];
				start_halfedges.push_back(he);
				start_halfedges.push_back(m.twinH[m.nextH[m.nextH[m.nextH[he]]]]);
			}

			else if (m.vertexDegree[vi] == 6 || m.vertexDegree[vi] == 8) { // singularity index = -1
				// expand towards all dominant edges (i.e. Emap[ei] = 0)
				int start_he = m.VH[vi];
				int he = start_he;
				do {
					if (m.Emap[m.HE[he]] == 0) {
						start_halfedges.push_front(he);
					}
					he = m.nextH[m.twinH[he]]; // attention, this only works if we are on a quad mesh
					if (he == -1) break;
					if(m.HV[he] != vi){cout << "Attention! HV[he] != vi! This cannot happen" << endl; throw; }
				} while (he != start_he);
			}

			else {
				cout << "Unknown vertex degree : " << m.vertexDegree[vi] << endl;
			}
		}
	}


	for (int he : start_halfedges)
		if (m.Emap[m.HE[he]] == 0)
			QuadMesh::expand_cutting_from_he(he, C, passed_vis, C.rows(), m.HE, m.HV, m.nextH, m.twinH, m.vertexDegree, m.boundaryVertices, stop_at_cut_intersections, true);
		else
			QuadMesh::expand_cutting_from_he(he, C, passed_vis, C.rows(), m.HE, m.HV, m.nextH, m.twinH, m.vertexDegree, m.boundaryVertices, true, true);



	///////////////////////////////////////////
	// --- create cuts from manually selected edges 
	const vector<int>& manually_selected_eis = np3dp->manually_selected_partition_eis;

	stop_at_cut_intersections = true;
	if (!manually_selected_eis.empty()) {
		//cout << "Also adding cuts to manually selected eis : ";
		for (int ei : manually_selected_eis) {
			//cout << ei << " , ";
			int he = m.EH(ei, 0);
			QuadMesh::expand_cutting_from_he(he, C, passed_vis, C.rows(), m.HE, m.HV, m.nextH, m.twinH, m.vertexDegree, m.boundaryVertices, stop_at_cut_intersections, true);
			he = m.nextH[m.twinH[m.nextH[m.twinH[he]]]];
			QuadMesh::expand_cutting_from_he(he, C, passed_vis, C.rows(), m.HE, m.HV, m.nextH, m.twinH, m.vertexDegree, m.boundaryVertices, stop_at_cut_intersections, true);
		}
		cout << endl;
	}

}


VectorXi Partitioning::create_Fmap(const VectorXi& D, const MatrixXi& F, const MatrixXi& EF, const MatrixXi& FE , const VectorXi& C)
{
	// --- color faces based on cuts (create Fmap)
	VectorXi Fmap;
	Fmap.setConstant(F.rows(), -1);

	int ID = -1;
	std::deque<int> Q;

	while (Fmap.minCoeff() == -1) { // while there exist uncolored faces

		if (Q.empty()) { // find one (any) uncolored face to expand from and store it in Q, also update current ID
			int fi = 0;
			while (Fmap[fi] != -1) {
				++fi;
			}
			Q.push_back(fi);
			ID = Fmap.maxCoeff() + 1;

			assert(Fmap[fi] == -1);
			Fmap[fi] = ID;
		}

		int fi = Q.front();
		Q.pop_front();

		for (int i = 0; i < D[fi]; ++i)
		{
			int ei = FE(fi, i);
			if (C[ei] == 0) { // if the ei is not marked as a cut
				int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
				if (nfi >= 0) { // if there exists a neighboring face
					if (Fmap[nfi] == -1) { // if neighbor face is uncolored
						Q.push_back(nfi);

						assert(Fmap[nfi] == -1);
						Fmap[nfi] = ID;
					}
				}
			}
		}
	}

	// --- Only accept partitions that have more than two faces
	/*int nP = Fmap.maxCoeff() + 1; // number of pieces
	VectorXi partitions_face_counts; partitions_face_counts.setZero(nP);
	for (int fi = 0; fi<F.rows(); ++fi)
		partitions_face_counts[Fmap[fi]] += 1;
	for (int ID=0; ID<nP; ++ID)
		if (partitions_face_counts[ID] < 2)
			for (int fi = 0; fi < F.rows(); ++fi)
				if (Fmap[fi] == ID)
					Fmap[fi] = -1;*/


	return Fmap;
}


namespace ConnectionsPrep
{
	void remove_pair(vector<pair<int,int>>& pairs, int id1, int id2)
	{
		pair<int, int> p1(id1, id2);
		pair<int, int> p2(id2, id1);
		int i = Helpers::element_index(pairs, p1);
		int j = Helpers::element_index(pairs, p2);
		if (i == -1 && j == -1) return;

		if (i >= 0)
			pairs.erase(pairs.begin() + i);
		else
			pairs.erase(pairs.begin() + j);
	}

	bool pair_in_vector(vector<pair<int,int>>& pairs, int id1, int id2)
	{
		pair<int, int> p1(id1, id2);
		pair<int, int> p2(id2, id1);
		return Helpers::in_vector(pairs, p1) || Helpers::in_vector(pairs, p2);
	}

	pair<int, int> get_pair(vector<pair<int,int>>& pairs, int id1, int id2)
	{
		pair<int, int> p1(id1, id2);
		if (Helpers::in_vector(pairs, p1))
			return p1;
		pair<int, int> p2(id2, id1);
		if (Helpers::in_vector(pairs, p2))
			return p2;
		cerr << "Attention! Could not find pair in pairs" << endl; throw;
	}

	// --- helper function
	bool fi_neighbors_attribute(int fi, int ID, const VectorXi& Attribute_map, const MatrixXi& EF, const MatrixXi& FE, const VectorXi& D, int& ei){
		// returns true if fi has a neighboring face nfi for which : Attribute_map[nfi] == ID
		// also stores in 'ei' the edge index through which the faces are neighboring
		for (int m = 0; m < D[fi]; ++m) {
			ei = FE(fi, m);
			const int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
			if (nfi >= 0)
				if (Attribute_map[nfi] == ID)
					return true;
		}
		return false;
	}

	void order_pieces_around_D8(vector<int>& pieces, const MatrixXi& F, const VectorXi& D, const VectorXi& Fmap, const MatrixXi& FE, const MatrixXi& EF)
	// considers the first piece as A, and orders the other 3 pieces around A
	{
		if (!pieces.size() == 4) { cout << "Wrong number of pieces while ordering around D8!" << endl; throw; }

		int p0 = pieces[0];
		int p1 = pieces[1];
		int p2 = pieces[2];
		int p3 = pieces[3];

		// find a direct neighbor of p0
		int neighbor = -1;
		for (int fi=0; fi<F.rows(); ++fi){
			if (Fmap[fi] == p0)
			{
				for (int k = 0; k < D[fi]; ++k) {
					int ei = FE(fi, k);
					int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
					if (Fmap[nfi] == p1 || Fmap[nfi] == p2 || Fmap[nfi] == p3){
						neighbor = Fmap[nfi];
						goto cnt;
					}
				}
			}
		}
		cnt:;
		if (neighbor == -1) { cout << "Could not find direct neighbor of p0" << endl; throw;}
	

		// find a direct neighbor of dn that is not p0
		int neighbor_next = 0;
		for (int fi=0; fi<F.rows(); ++fi){
			if (Fmap[fi] == neighbor)
			{
				for (int k = 0; k < D[fi]; ++k) {
					int ei = FE(fi, k);
					int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
					if (Fmap[nfi] != p0 && Fmap[nfi] != neighbor){
						if (Fmap[nfi] == p1 || Fmap[nfi] == p2 || Fmap[nfi] == p3){
							neighbor_next = Fmap[nfi];
							goto cnt1;
						}
					}
				}
			}
		}
		cnt1:;
		if (neighbor_next == -1) { cout << "Could not find direct neighbor of dn" << endl; throw;}


		// get last remaining piece
		int lp = -1;
		if (p1 != neighbor && p1 != neighbor_next)      lp = p1;
		else if (p2 != neighbor && p2 != neighbor_next) lp = p2;
		else											lp = p3;


		// save found order
		pieces[1] = neighbor;
		pieces[2] = neighbor_next;
		pieces[3] = lp;
	}
}


void Partitioning::separate_pieces(Np3dpContext* np3dp, const VectorXi& Emap_original, const MatrixXd& UV_original)
{
	VectorXi& Fmap = np3dp->Fmap;
	const VectorXi& ERibMap_original = np3dp->quad_mesh.ERibMap;

	if (Emap_original.rows() != ERibMap_original.rows()) { cerr << "Attention! np3dp->Emap.rows() != np3dp->ERibMap.rows()" << endl; throw; }

	int nP = Fmap.maxCoeff() + 1; // number of pieces

	np3dp->pieces.clear();
	np3dp->pieces.resize(nP);

	cout << endl;
	for (int ID=0; ID<nP; ++ID) { // ID: piece id
		cout << "Separating quad mesh piece with ID = " << ID << endl;

		MatrixXd& V = np3dp->pieces[ID].V; // current piece vertices
		MatrixXi& F = np3dp->pieces[ID].F; // current piece faces
		VectorXi& D = np3dp->pieces[ID].D; // current piece face degrees
		VectorXi& V_to_Vquad_original = np3dp->pieces[ID].V_to_Vquad_original; // keep track of correspondence of current piece's vertex to original quad mesh
		VectorXi& Fs_to_Fquad_original = np3dp->pieces[ID].F_to_Fquad_original;
		VectorXi& originalBoundaryFaces = np3dp->pieces[ID].originalBndryFaces;

		VectorXi V_hash; V_hash.setConstant(np3dp->quad_mesh.V.rows(), -1); // -1 if vertex is unvisited, index of vertex in new V matrix if visited
		VectorXi F_visited; F_visited.setZero(np3dp->quad_mesh.F.rows()); // keep track of which faces have already been visited

		// --- find starting fi (any face that belongs to the current piece)
		int start_fi = -1;
		for (int fi = 0; fi < np3dp->quad_mesh.F.rows(); ++fi) {
			if (Fmap[fi] == ID) {
				start_fi = fi; break;
			}
		}
		if (start_fi == -1) {cerr<< "Could not find a starting face!" << endl; throw;} // make sure a start_fi was found

		// --- expand from start_fi to visit the whole piece
		std::deque<int> Q;
		Q.push_back(start_fi);
		F_visited[start_fi] = 1;

		while(!Q.empty()) {
			int fi = Q.front();
			Q.pop_front();

			// add new face to F matrix, and new vertices to V matrix
			D.conservativeResize(F.rows() + 1);
			Fs_to_Fquad_original.conservativeResize(F.rows() + 1);
			originalBoundaryFaces.conservativeResize(F.rows() + 1);
			F.conservativeResize(F.rows() + 1, np3dp->quad_mesh.F.cols()); // add new face

			const int fc = F.rows() - 1;

			D(fc) = np3dp->quad_mesh.D[fi];


			Fs_to_Fquad_original(fc) = fi;
			originalBoundaryFaces(fc) = np3dp->quad_mesh.boundaryFaces[fi];


			for (int m=0; m<np3dp->quad_mesh.D[fi]; ++m) {
				int vi = np3dp->quad_mesh.F(fi, m);

				if (V_hash[vi] == -1) { // if vertex is unvisited, set hash to index in new V array. Otherwise that vertex + hash are already set
					V_to_Vquad_original.conservativeResize(V.rows() + 1);
					V.conservativeResize(V.rows() + 1, 3); // add new vertex

					V.row(V.rows() - 1) = np3dp->quad_mesh.V.row(vi);
					V_to_Vquad_original[V.rows() - 1] = vi;
					V_hash[vi] = V.rows() - 1;
				}

				F(fc, m) = V_hash[vi];

				// --- append eligible face neighbors in Q
				int ei = np3dp->quad_mesh.FE(fi, m);
				int nfi = np3dp->quad_mesh.EF(ei, 0) == fi ? np3dp->quad_mesh.EF(ei, 1) : np3dp->quad_mesh.EF(ei, 0);

				if (nfi >= 0) {
					if ((Fmap[nfi] == ID) && F_visited[nfi] == 0) {
						Q.push_back(nfi);
						F_visited[nfi] = 1;
					}
				}
			}
		}


		// TODO: remove --->
		// cout << "-------------- ATTENTION! --------------" << endl;
		// cout << "-------------- Flipping Ftype! --------------" << endl;
		// Ftype *= -1;
		// cout << "-------------- ATTENTION! --------------" << endl;
		// TODO: <--- remove

		// --- get uv coords from original mesh
		MatrixXd& UV = np3dp->pieces[ID].UV;
		UV.setZero(V.rows(), 2);
		for (int vi = 0; vi < V.rows(); ++vi)
			UV.row(vi) = UV_original.row(V_to_Vquad_original[vi]);

		// --- get FixedVis from original mesh (we are not really using this value though)
		// VectorXi FixedVis; FixedVis.setZero(V.rows());
		// todo: use V_to_Vquad_original to get true values (only if they are actually needed)

		// --- update Emap and Erib
		MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
		hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
		VectorXi& Emap = np3dp->pieces[ID].Emap;
		VectorXi& Erib = np3dp->pieces[ID].ERibMap;
		get_Emap_and_Erib_from_original_quad_mesh(Emap, Erib, Emap_original, ERibMap_original, V_to_Vquad_original, EF, EV, np3dp->quad_mesh.V, np3dp->quad_mesh.EV);

		// TODO: cleanup individual piece?
	}
}


void Partitioning::get_Emap_and_Erib_from_original_quad_mesh(VectorXi& Emap, VectorXi& Erib, const VectorXi& Emap_original, const VectorXi& ERibMap_original, 
	const VectorXi& V_to_Vquad_original, const MatrixXi& EF, const MatrixXi& EV, const MatrixXd& Vquad_original, const MatrixXi& EVquad_original)
{
	// --- get partitioning VE_original from original (i.e. unpartitioned) mesh's edges
	vector<vector<int>> VE_original(Vquad_original.rows()); // VE map of original quad mesh
	for (int ei = 0; ei < EVquad_original.rows(); ++ei) {
		int v0 = EVquad_original(ei, 0);
		int v1 = EVquad_original(ei, 1);
		VE_original[v0].push_back(ei);
		VE_original[v1].push_back(ei);
	}

	// --- update Emap and Erib
	Emap.setConstant(EF.rows(), -1);
	Erib.setZero(EF.rows());
	for (int ei = 0; ei < EF.rows(); ++ei) {
		int first_v = EV(ei, 0);
		int second_v = EV(ei, 1);
		int v0 = V_to_Vquad_original[EV(ei, 0)]; // vertex indices of original mesh
		int v1 = V_to_Vquad_original[EV(ei, 1)];

		if (v0 == v1) // this should only happen if this edge doesn't exist at all in the original mesh (can only happen after the cutting of faces for the removal of gaps)
		{
			Emap[ei] = -1;
			Erib[ei] = 0;
		}
		else
		{
			// find common element in two edge lists (i.e. common edge between v0 and v1)
			vector<int> common_elements(VE_original[v0].size());// + VE_original[v1].size());
			vector<int>::iterator it, end;
			end = set_intersection(VE_original[v0].begin(), VE_original[v0].end(),
				VE_original[v1].begin(), VE_original[v1].end(),
				common_elements.begin());

			assert(*end == 0); // make sure there's only 1 common element
			int common_edge = common_elements[*end];
			Emap[ei] = Emap_original[common_edge];
			Erib[ei] = ERibMap_original[common_edge];
		}
	}
}



namespace Partitioning_Helpers
{
	double get_sum_of_angles_from_strip(const vector<int>& strip_fis, const MatrixXd& GuidingDirs) // summing up angles in the green direction
	{
		double angle_sum = 0;
		for (int i=0; i<strip_fis.size()-1; ++i)
		{
			int fi  = strip_fis[i];
			int nfi = strip_fis[i+1];
			const RowVector3d& green_dir1 = GuidingDirs.row(fi);
			const RowVector3d& green_dir2 = GuidingDirs.row(nfi);
			angle_sum += acos(green_dir1.dot(green_dir2)); // for x = [-1, 1] acos returns [0, pi] radians
		}
		return angle_sum;
	}


	void get_piece_max_angle_variation_in_subdominant_direction(const Piece& piece, const MatrixXi& strips, const VectorXi& strips_directions, const vector<vector<int>>& fis,
		double& blue_strips_max_angle_sum, int& bsi, double& green_strips_max_angle_sum, int& gsi)
	{
		// check each vector with each other vector, and keep the max angle 
		// for each blue strip, examine green direction
		// for each green strip, examine (ALSO) green direction
		// bf1, bf2 and  gf1, gf2 are the faces where the min dot was found. bsi and gsi are the strip indices with the min dot

		blue_strips_max_angle_sum = 0;
		green_strips_max_angle_sum = 0;


		for (int si=0; si<strips.rows(); ++si)
		{
			double angle_sum = get_sum_of_angles_from_strip(fis[si], piece.GuidingDirs);

			if (strips_directions[si] == 0) // blue strip
			{
				if (angle_sum > blue_strips_max_angle_sum){
					blue_strips_max_angle_sum = angle_sum;
					bsi = si;
				}
			}
			else // green strip
			{
				if (angle_sum > green_strips_max_angle_sum){
					green_strips_max_angle_sum = angle_sum;
					gsi = si;
				}
			}
		}
	}


	int find_longest_strip_in_direction(int dir, const VectorXi& strips_directions, const VectorXi& strip_sums)
	{
		int longest_si = -1; int max_size = 0; // the strip that has most faces
		for (int si = 0; si < strips_directions.rows(); ++si)
		{
			if (strips_directions[si] == dir)
			{
				if (max_size < strip_sums[si]) {
					longest_si = si;
					max_size = strip_sums[si];
				}
			}
		}
		return longest_si;
	}
}


void Partitioning::identify_snake_strips(Np3dpContext* np3dp, VectorXd& angle_sums) // works on the unpartitioned quad mesh
{
	QuadMesh& q = np3dp->quad_mesh;

	if (q.EV.rows() != q.Emap.rows()) { throw invalid_argument("The Emap of the piece has the wrong number of rows"); }

	// --- get all strips (only in the dominant direction)
	vector<array<vector<int>, 2>> strips_vis;
	vector<bool> strips_are_closed;

	Strips::get_all_strips(q.V, q.F, q.D, q.boundaryFaces, q.Emap, q.strip_networks.strips, q.strip_networks.StoD, strips_vis, strips_are_closed, q.strip_networks.StoEis, 
		q.strip_networks.StoFis, q.strip_networks.FtoS, q.strip_networks.VtoS, {0});

	// debug -->
	// Helpers::write_matrix_to_txt(np3dp->strips, DATA_PATH + np3dp->output_folder + "piece_strips.txt");
	// Helpers::write_vector_to_txt(np3dp->StoD, DATA_PATH + np3dp->output_folder + "piece_strips_directions.txt");
	// <---

	auto get_face_dir = [](int fi, const MatrixXd& V, const VectorXi& D, const MatrixXi& FE, const MatrixXi& EV, const VectorXi& Emap, int emap_value)->RowVector3d
	{
		RowVector3d dir(0, 0, 0);
		for (int i=0; i<D[fi]; ++i)
		{
			int ei = FE(fi, i);
			if (Emap[ei] == emap_value)
			{
				Vector3d current_dir = V.row(EV(ei, 0)) - V.row(EV(ei, 1));
				if (current_dir.dot(dir) < 0)
					current_dir *= -1;
				dir += current_dir;
			}
		}
		dir.normalize();
		return dir;
	};


	angle_sums.setZero(q.strip_networks.strips.rows());
	for (int si=0; si<q.strip_networks.strips.rows(); ++si)
	{
		MatrixXd green_dirs(q.strip_networks.StoFis[si].size(), 3);
		if (q.strip_networks.StoD[si] == 0)
		{
			for (int i=0; i<q.strip_networks.StoFis[si].size() - 1; ++i)
			{
				int fi  = q.strip_networks.StoFis[si][i];
				int nfi = q.strip_networks.StoFis[si][i+1];

				RowVector3d green_dir1;
				if (i == 0){
					green_dir1 = get_face_dir(fi,  q.V, q.D, q.FE, q.EV, q.Emap, 0);
					green_dirs.row(i) = green_dir1;
				}
				else{
					green_dir1 = green_dirs.row(i);
				}
				
				RowVector3d green_dir2 = get_face_dir(nfi, q.V, q.D, q.FE, q.EV, q.Emap, 0);
				
				if (green_dir1.dot(green_dir2) < 0) // make sure direction is consistent (assuming relatively 'smooth' mesh so that the transition from one face to the next has ange less than PI)
					green_dir2 *= -1;
				green_dirs.row(i+1) = green_dir2;

				if ((green_dir1 - green_dir2).squaredNorm() > 1e-8) // if two vectors are identical the acos(1) returns nan
					angle_sums[si] += acos(green_dir1.dot(green_dir2)); // for x = [-1, 1] acos returns [0, pi] radians
			}
		}
		// Helpers::write_matrix_to_txt(green_dirs, DATA_PATH + np3dp->output_folder + "green_dirs.txt");
	}
}