//
//  Created by Ioanna Mitropoulou.
//

#include "piece.h"
#include <directional/polygonal_edge_topology.h>
#include "quad_mesh.h"
#include <hedra/triangulate_mesh.h>

void Piece::get_avg_directions_and_extents(Vector3d& extents, Matrix3d& directions) const
{
    if(F.rows() != GuidingDirs.rows() || F.rows() != faceNormals.rows()) throw invalid_argument("The GuidingDirs or per-face Normals do not have the same dimensions as F");

    // --- get directions
    RowVector3d up_dir(0,0,0), forward_dir(0,0,0);
    for (int fi=0; fi< F.rows(); ++fi)
    {
        up_dir += GuidingDirs.row(fi);
        RowVector3d n = faceNormals.row(fi); //block(fi, 0, 1, 3);
        RowVector3d d = GuidingDirs.row(fi);
        forward_dir += n.cross(d);
    }
    up_dir.normalize();
    forward_dir.normalize();

    /* Ioanna note: I think the following is a more correct way for getting the forward_dir 
     * 				// get green direction (average of two green edges that points in the direction of forward dir)
					RowVector3d blue_dir = piece.GuidingDirs.row(fi);
					RowVector3d n = piece.N.row(fi);
					RowVector3d forward_dir = n.cross(blue_dir);

					RowVector3d green_dir(0, 0, 0);
					for (int k = 0; k < piece.D[fi]; ++k){
						int ei = FE(fi, k);
						if (piece.Emap[ei] == 1){
							RowVector3d dir = piece.V.row(EV(ei, 1) - EV(ei, 0)).normalized();
							if (dir.dot(forward_dir) < 0)
								dir *= -1;
							green_dir += dir;
						}
					}
					green_dir.normalize();
     */

    // set each axis at each column of direction
    directions.col(0) = up_dir;
    directions.col(1) = forward_dir;
    directions.col(2) = up_dir.cross(forward_dir);
    directions.transposeInPlace(); // since it's a rotation matrix, the transpose is also the inverse of the transformation

    // --- get extents
    MatrixXd V_in_directions_coordinate_system = (directions * V.transpose()).transpose(); // Rotate points from world to the coordinate system of the eigenvectors // Pcopy is transposed (i.e. 2xnP)
    //Helpers::write_matrix_to_txt(V_in_directions_coordinate_system, "C:/dev/np3dp_v2/data/curvy_surface_2/output/V_in_directions_coordinate_system.txt");

	extents.setZero(3);
    for (int i = 0; i < 3; ++i)
    {
        extents[i] = V_in_directions_coordinate_system.col(i).maxCoeff() - V_in_directions_coordinate_system.col(i).minCoeff();
    }
}


MatrixXd Piece::get_per_face_normals(const MatrixXd& N_ref)
{
    MatrixXi T; VectorXi TF;
    hedra::triangulate_mesh(D, F, T, TF);
    MatrixXd N = QuadMesh::per_face_normals(V, T, F, TF);
    int f0 = 0;
    while (F_to_Fquad_original[f0] == -1)
        ++f0;
    int f0_original =F_to_Fquad_original[f0];
    if (N_ref.row(f0_original).dot(N.row(f0)) < 0)
        N = N * (-1);

    return N;
}

void Piece::serialize(int piece_id, const std::string& serialize_folder) const 
{
    const std::string filename_start = DATA_PATH + serialize_folder + "Piece_" + std::to_string(piece_id);

    igl::serialize(V,                   filename_start+ "_V.igl");
    igl::serialize(F,                   filename_start + "_F.igl");
    igl::serialize(D,                   filename_start + "_D.igl");
    igl::serialize(faceNormals,         filename_start + "_N.igl");
    igl::serialize(Emap,                filename_start + "_Emap.igl");
    igl::serialize(ERibMap,             filename_start + "_ERibMap.igl");
    igl::serialize(V_to_Vquad_original, filename_start + "_V_to_Vquad_original.igl");
    igl::serialize(F_to_Fquad_original, filename_start + "_F_to_Fquad_original.igl");
	igl::serialize(originalBndryFaces,  filename_start + "_originalBoundaryFacesQuad.igl");
    igl::serialize(UV,                  filename_start + "_UV.igl");
    igl::serialize(GuidingDirs,         filename_start + "_GuidingDirs.igl");
    igl::serialize(nP_per_strip,         filename_start + "_nP_per_strip.igl");

    // currently PathCollection is NOT serialized
}

void Piece::deserialize(int piece_id, const std::string& serialize_folder)
{
    const std::string filename_start = DATA_PATH + serialize_folder + "Piece_" + std::to_string(piece_id);

    igl::deserialize(V,                   filename_start+ "_V.igl");
    igl::deserialize(F,                   filename_start + "_F.igl");
    igl::deserialize(D,                   filename_start + "_D.igl");
    igl::deserialize(faceNormals,         filename_start + "_N.igl");
    igl::deserialize(Emap,                filename_start + "_Emap.igl");
    igl::deserialize(ERibMap,             filename_start + "_ERibMap.igl");
    igl::deserialize(V_to_Vquad_original, filename_start + "_V_to_Vquad_original.igl");
    igl::deserialize(F_to_Fquad_original, filename_start + "_F_to_Fquad_original.igl");
	igl::deserialize(originalBndryFaces,  filename_start + "_originalBoundaryFacesQuad.igl");
    igl::deserialize(UV,                  filename_start + "_UV.igl");
    igl::deserialize(GuidingDirs,         filename_start + "_GuidingDirs.igl");
    igl::deserialize(nP_per_strip,         filename_start + "_nP_per_strip.igl");
}


namespace GuidingDirsHelpers
{
    RowVector3d get_undirected_face_guiding_dir(int fi, const MatrixXd& V, const VectorXi& D, const MatrixXi& FH, const VectorXi& HE, const VectorXi& HV, const VectorXi& nextH,
        const VectorXi& Emap, const MatrixXd& faceNormals) {
        RowVector3d guiding_dir(0.0, 0.0, 0.0);
        int count_found_edges = 0;
        for (int m = 0; m < D[fi]; ++m) {
            int he = FH(fi, m);
            int ei = HE[he];
            int v1 = HV[he];
            int v2 = HV[nextH[he]];
            if (Emap[ei] == 0) { // dominant edge (uv_diff = 0)
                int dir = count_found_edges == 0 ? 1 : -1;
                guiding_dir += (V.row(v2) - V.row(v1)).normalized() * dir;
                ++count_found_edges;
            }
        }
        if (count_found_edges == 0) { cerr << "While creating guiding directions, could not find any edges with Emap[ei]==0 (i.e. dominant) on fi = " << fi << endl; return RowVector3d(0,0,0);}
        guiding_dir = guiding_dir.cross(RowVector3d(faceNormals(fi, 0), faceNormals(fi, 1), faceNormals(fi, 2))); // rotate 90 degrees
        guiding_dir.normalize();
        double n = guiding_dir.squaredNorm();
        if (!(0.9 < n && n < 1.1)) { cerr << "Guiding direction for fi = " << fi << " has norm != 1. This should never happen." << endl; throw; }
        return guiding_dir.normalized();
    }
}


MatrixXd Piece::get_guiding_directions()
{
    MatrixXd guiding_dirs;

	// edge topology data
    MatrixXi EV, FE, EF, EFi;
	MatrixXd FEs;
    VectorXi innerEdges;

	// halfedge data
    MatrixXi EH, FH;
	VectorXi VH, HE, HF, HV, nextH, prevH, twinH;

	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);


    // --- get "virtual" cuts from any wedge singularities on the interior of the piece
	VectorXi vertexDegree = QuadMesh::mesh_vertex_degrees(V, F, D);
    VectorXi boundaryVertices = QuadMesh::vertices_on_boundaries(V, EF, EV);
    VectorXi boundaryFaces = QuadMesh::faces_on_boundaries(F, EF);


	hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);


    VectorXi VirtualCuts;
    VirtualCuts.setZero(EV.rows()); // cuts (#E x 1) 0: not a cut, 1: cut
    deque<int> start_halfedges;
    for (int vi = 0; vi < V.rows(); ++vi)
    {
        if (vertexDegree[vi] != 4 && boundaryVertices[vi] != 1) {
            cout << "Starting virtual cut from vi = " << vi << " , with vertex degree : " << vertexDegree[vi] << endl;

            if (vertexDegree[vi] == 2) { // singularity index = 1
                int he = VH[vi];
                if (Emap[HE[he]] == 0) {
                    start_halfedges.push_front(he);
                    VirtualCuts[HE[twinH[nextH[nextH[nextH[he]]]]]] = 1;
                }
                else {
                    start_halfedges.push_front(twinH[nextH[nextH[nextH[he]]]]);
                    VirtualCuts[HE[he]] = 1;
                }
            }
        }
    }
    for (int he : start_halfedges)
        QuadMesh::expand_cutting_from_he(he, VirtualCuts, vector<int>(), VirtualCuts.rows(), HE, HV, nextH, twinH, vertexDegree, boundaryVertices, false, true);
    //cout << "VirtualCuts.sum() = " << VirtualCuts.sum() << endl;


    // --- set guiding directions
    int f_start = 0; // arbitrary selection of first face
    while (D[f_start] != 4 && f_start < F.rows() - 1)
        ++f_start;

    guiding_dirs.setZero(F.rows(), 3);
    guiding_dirs.row(f_start) = GuidingDirsHelpers::get_undirected_face_guiding_dir(f_start, V, D, FH, HE, HV, nextH, Emap, faceNormals); // get first forward direction (the sign doesn't matter, it will be kept consistent for all quads visited after this one)

    std::deque<int> P = { f_start };
    VectorXi visited; visited.setZero(F.rows()); visited[f_start] = 1;

    while (!P.empty()) {
        int fi = P.front();
        P.pop_front();
        double n = guiding_dirs.row(fi).squaredNorm(); assert(0.9 < n && n < 1.1); // make sure the guiding dir has been set and normalized

        for (int m = 0; m < D[fi]; ++m) {
            int ei = FE(fi, m);
            if (VirtualCuts[ei] == 1) continue; // do not cross any virtual cut
            int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
            if (nfi >= 0) {
                if (!visited[nfi])
                {
                    // set nfi guiding direction based on fi
                    RowVector3d new_guiding_dir = GuidingDirsHelpers::get_undirected_face_guiding_dir(nfi, V, D, FH, HE, HV, nextH, Emap, faceNormals);
                    double n = new_guiding_dir.squaredNorm();
                    if (n < 0.1) // if we couldnt find a guiding dir, then just use the neighbor's
                        new_guiding_dir = guiding_dirs.row(fi);
                    double p = guiding_dirs.row(fi).squaredNorm(); assert(0.9 < p && p < 1.1);

                	if (abs(guiding_dirs.row(fi).dot(new_guiding_dir)) < 0.2) // if we are on a very sharp corner, then do not consider this neighboring. Instead continue to find another route 
                        continue; // on very sharp corners undesirable sign flips can occur

                    if (guiding_dirs.row(fi).dot(new_guiding_dir) < 0) // if they are pointing in opposite directions
                        new_guiding_dir *= -1;
                    guiding_dirs.row(nfi) = new_guiding_dir;

                    // add nfi to queue
                    P.push_back(nfi);
                    visited[nfi] = 1;
                }
            }
        }
    }

    return guiding_dirs;
}


int Piece::piece_ei_to_original_mesh_ei(int piece_v0, int piece_v1, QuadMesh& quad_mesh) const
{
    const int v0 = V_to_Vquad_original[piece_v0];
    const int v1 = V_to_Vquad_original[piece_v1];

    for (int ei = 0; ei < quad_mesh.EV.rows(); ++ei) {
        if (quad_mesh.EV(ei, 0) == v0 && quad_mesh.EV(ei, 1) == v1 || 
            quad_mesh.EV(ei, 0) == v1 && quad_mesh.EV(ei, 1) == v0 )
        return ei;
    }
    return -1;
}