//
// Created by Ioanna Mitropoulou on 07.01.22.
//

#include "quad_mesh.h"

#include <directional/mesh_function_isolines.h>
#include <directional/setup_mesh_function_isolines.h>
#include <hedra/triangulate_mesh.h>
#include <directional/polygonal_edge_topology.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/matlab/matlabinterface.h>
#include <hedra/planarity.h>

#include "partitioning.h"
#include "strips.h"


void QuadMesh::generate_from_integration(const directional::TriMesh& mesh, const directional::TriMesh& meshCut, const directional::IntegrationData& intData, const MatrixXd& triangle_mesh_UVcoords)
{ 
	// --- setup
	cout << "Setting up quad mesh generation" << endl;
	directional::MeshFunctionIsolinesData mfiData;
	directional::setup_mesh_function_isolines(meshCut, intData, mfiData);

	// --- meshing
	directional::mesh_function_isolines(mesh, mfiData, intData.verbose, V, D, F);
	cout << "Generated quad mesh with " << F.rows() << " quads." << endl;

	// --- turn pairs of triangles on singularities into quads
	hedra::triangulate_mesh(D, F, T, TF);
	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);
	add_extra_vertex_on_triangle_faces();

	// get UVcoords
	get_UVcoords_on_quad_vertices(meshCut, triangle_mesh_UVcoords);

	update_quad_mesh_data();

	cleanup();

    cout << "Creating Emap for quad mesh" << endl;
    create_quad_Emap(false);

    // saving mesh data
    save_data();
}

void QuadMesh::update_quad_mesh_data() {

    hedra::triangulate_mesh(D, F, T, TF);
	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

    faceNormals = per_face_normals(V, T, F, TF);
    vertexDegree = mesh_vertex_degrees(V, F, D);
    boundaryVertices = vertices_on_boundaries(V, EF, EV);
    boundaryFaces = faces_on_boundaries(F, EF);
}

double QuadMesh::get_avg_len_edge_type(int edge_type) const
{
    int count = 0;
    double d = 0;
	for (int ei=0; ei<Emap.size(); ++ei) {
        if (Emap[ei] == edge_type) {
	        int v0 = EV(ei, 0);
            int v1 = EV(ei, 1);
            d += (V.row(v0) - V.row(v1)).norm();
            ++count;
        }
	}
    return d / static_cast<double>(count);
}

VectorXd QuadMesh::face_areas(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const VectorXi& TF)
{
	VectorXd At;
    igl::doublearea(V, T, At);
    At *= 0.5;

    VectorXd A;
	A.setZero(F.rows());
    for (int ti=0; ti<T.rows(); ++ti) {
        int fi = TF[ti];
        A[fi] += At[ti];
    }
    return A;
}

bool QuadMesh::has_updated_mesh_data() const
{
    if (V.rows() != UVcoords.rows()) { cout << "Vquad.rows() != UVcoordsQuad.rows()" << endl; return false; }
    if (Emap.rows() != EV.rows()) { cout << "Emap.rows() != EV.rows()" << endl; return false; }
    if (V.rows() != boundaryVertices.rows()) { cout << "Vquad.rows() != boundaryVerticesQuad.rows()" << endl; return false; }
    if (F.rows() != boundaryFaces.rows()) { cout << "Fquad.rows() != boundaryFacesQuad.rows()" << endl; return false; }
    return true;
}

void QuadMesh::remember_previous_data()
{
	V_prev = V;
	F_prev = F;
	D_prev = D;
	UV_prev = UVcoords;
	Emap_prev = Emap;
	has_previous_data = true;
}

void QuadMesh::restore_previous_data()
{
    // set prev data
    V = V_prev;
    F = F_prev;
    D = D_prev;
    UVcoords = UV_prev;
    Emap = Emap_prev;

    invalidate_previous_data();

    // --- update np3dp data
    update_quad_mesh_data();
    if (!has_updated_mesh_data()) { cout << "QuadMeshEditing::undo(Np3dp& np3dp) finished without updated mesh data" << endl; throw; }

    save_data();
}

void QuadMesh::invalidate_previous_data()
{
	// delete prev data
    has_previous_data = false;
    V_prev.setZero(0, 0);
    F_prev.setZero(0, 0);
    D_prev.setZero(0);
    UV_prev.setZero(0, 0);
    Emap_prev.setZero(0);
}


VectorXi QuadMesh::mesh_vertex_degrees(const MatrixXd& V, const MatrixXi& F, const VectorXi& D)
{
    VectorXi vertexDegree;
    int nV = V.rows();
	vertexDegree.setZero(nV);
	for (int fi = 0; fi < F.rows(); ++fi) {
        for (int m = 0; m < D[fi]; ++m) {
            int vi = F(fi, m);
            if (vi >= 0) {
	            vertexDegree[vi] += 1;
            }
            else {
                cout << "Attention! On fi = " << fi << " with D[fi] = " << D[fi] << ", there is a vi = -1 in position m = " << m << ". F.row(fi) = " << F.row(fi) << endl; throw;
            }
        }
    }
    return vertexDegree;
}


VectorXi QuadMesh::vertices_on_boundaries(const MatrixXd& V, const MatrixXi& EF, const MatrixXi& EV) {
    VectorXi BoundaryVis;
	int nV = V.rows();
    BoundaryVis.setZero(nV);

    for (int ei=0; ei<EF.rows(); ++ei)
    {
	    int f1 = EF(ei, 0);
        int f2 = EF(ei, 1);
        if (f1 == -1 || f2 == -1) // then edge is on boundary
        {
            int v1 = EV(ei, 0);
            int v2 = EV(ei, 1);
            BoundaryVis[v1] = 1;
            BoundaryVis[v2] = 1;
        }
    }
    return BoundaryVis;
}


VectorXi QuadMesh::faces_on_boundaries(const MatrixXi& F, const MatrixXi& EF)
{
    VectorXi BoundaryFis;
    BoundaryFis.setZero(F.rows());

    for (int ei = 0; ei < EF.rows(); ++ei)
    {
        int f1 = EF(ei, 0);
        int f2 = EF(ei, 1);
        if (f1 == -1 || f2 == -1) // then edge is on boundary
        {
            int bfi = f1 == -1 ? f2 : f1;
            BoundaryFis[bfi] = 1;
        }
    }
    return BoundaryFis;
}


bool QuadMesh::is_boundary_edge(int ei, const MatrixXi& EH){
    return EH(ei, 0) == -1 || EH(ei, 1) == -1;
}


bool QuadMesh::is_boundary_face(int f0, const VectorXi& twinH, const VectorXi& D, const MatrixXi& FH){
    for (int m = 0; m < D[f0]; ++m) {
        int he = FH(f0, m);
        if (twinH[he] == -1)
            return true;
    }
    return false;
}


MatrixXd QuadMesh::get_lines_matrix_from_cuts(const VectorXi& C) const
{
    int n = C.sum();
    MatrixXd lines(n, 6);
    int i = 0;
    for (int ei = 0; ei < EV.rows(); ++ei) {
        if (C[ei]) {
            lines.block(i, 0, 1, 3) = V.row(EV(ei, 0));
            lines.block(i, 3, 1, 3) = V.row(EV(ei, 1));
            ++i;
        }
    }
    return lines;
}


void QuadMesh::create_quad_Emap(bool remember_existing_emap_values)
{
    assert(FE.rows() == F.rows());
    int nE = EF.rows();

    // --- Remember Emap for one ei (that is not on the boundary)
    int e_remember_v1=-1, e_remember_v2=-1, e_remember_Emap=-1;
    if (remember_existing_emap_values && Emap.rows() == nE)
    {
        int ei = 0;
        while (true) {
            if (Emap[ei] != -1)
                break;
            ++ei;
        }
        e_remember_v1 = EV(ei, 0);
        e_remember_v2 = EV(ei, 1);
        e_remember_Emap = Emap[ei];
    }

    Emap.setConstant(nE, -1);

    // traverse all faces and 'color' edges as either 0 or 1
    int nF = F.rows();
    VectorXi visited;
    visited.setZero(nF);

    int f0 = 0; // arbitrary selection of fi=0 as starting face
    while (D[f0] != 4 || is_boundary_face(f0, twinH, D, FH)) { // Make sure that the first selected face is a quad
        ++f0;
    }
    //cout << "Starting face for Emap: f0 = " << f0 << endl;

    std::deque<int> Q = { f0 };
    visited[f0] = 1;

    // --- starting face : arbitrary coloring of its edges
    Emap[FE(f0, 0)] = 0;
    Emap[FE(f0, 1)] = 1;
    Emap[FE(f0, 2)] = 0;
    Emap[FE(f0, 3)] = 1;

    // --- expand and color other faces based on that
    while (!Q.empty())
    {
        int fi = Q.front();
        Q.pop_front();

        for (int i = 0; i < 4; ++i)
        {
            int ei = FE(fi, i);
            int nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
            if (nfi > 0) {
                if (D[nfi] == 4) { // only for quads
                    if (!visited[nfi]) {

                        // find which of the 4 edges of nfi is already colored
                        int colored_edge = -1;
                        int color = -1;
                        for (int k = 0; k < 4; ++k) {
                            int ei = FE(nfi, k);
                            if (Emap[ei] != -1) {
                                colored_edge = k;
                                color = Emap[ei];
                            }
                        }
                        if (colored_edge == -1) 
                            throw invalid_argument("Visited face without ANY colored edge!"); // there should always be at least one colored edge
                        

                        // color all edges of nfi
                        for (int k = 1; k < 4; ++k) {
                            int ei = FE(nfi, (colored_edge + k) % 4);
                            int new_color = k % 2 == 0 ? color : static_cast<int>(!static_cast<bool>(color));
                            if (Emap[ei] != -1)
                            {
                                if (Emap[ei] != new_color) {
                                    cerr << "Attention, color conflict at edge ei = " << ei << " : color = " << color << " , new_color = " << new_color << " ,  vi1 = " << EV(ei, 0) << " , vi2 = " << EV(ei, 1) << endl;
                                    cerr << "v1 coords = " << V.row(EV(ei, 0)) << endl;
                                    cerr << "v2 coords = " << V.row(EV(ei, 1)) << endl;
                                    //throw; // TODO: we should never get here, but it can happen after rewiring (in the areas where additional edges need to be added next to the boundaries)
                                }
                            }
                            else {
                                Emap[ei] = new_color;
                            }
                        }

                        // add nfi to Q
                        Q.push_back(nfi);
                        visited[nfi] = 1;
                    }
                }
            }
        }
    }


    auto color_from_nextHs = [](int ei, const MatrixXi& EH, const VectorXi& twinH, const VectorXi& nextH, const VectorXi& HE, const VectorXi& Emap)->int {
        for (int he : {EH(ei, 0), EH(ei, 1)}) {
            if (he >= 0) {
                he = twinH[he]; if (he == -1) continue;
                he = nextH[he]; if (he == -1) continue;
                he = twinH[he]; if (he == -1) continue;
                he = nextH[he]; if (he == -1) continue;
                int other_ei = HE[he]; // if you find extension
                if (Emap[other_ei] != -1) // if the extension has a set Emap
                    return Emap[other_ei];
            }
        }
        return -1;
    };

    auto color_from_prevHs = [](int ei, const MatrixXi& EH, const VectorXi& twinH, const VectorXi& prevH, const VectorXi& HE, const VectorXi& Emap)->int {
        for (int he : {EH(ei, 0), EH(ei, 1)}) {
            if (he >= 0) {
                he = prevH[he]; if (he == -1) continue;
                he = twinH[he]; if (he == -1) continue;
                he = prevH[he]; if (he == -1) continue;
                int other_ei = HE[he]; // if you find extension
                if (Emap[other_ei] != -1) // if the extension has a set Emap
                    return Emap[other_ei];
            }
        }
        return -1;
    };

    auto color_from_neighboring_edges = [](int ei, const MatrixXi& EH, const VectorXi& prevH, const VectorXi& nextH, const VectorXi& HE, const VectorXi& Emap)->int {
        int he0 = EH(ei, 0);
        int he1 = EH(ei, 1);
        if (he0 >= 0 && he1 >= 0) {

            // get color from f1
            int c0 = -1;
            int next_e0 = HE[nextH[he0]];
            int prev_e0 = HE[prevH[he0]];
            if (Emap[next_e0] != -1)
                c0 = int(!bool(Emap[next_e0]));
            else if (Emap[prev_e0] != -1)
                c0 = int(!bool(Emap[prev_e0]));

            // get color from f2
            int c1 = -1;
            int next_e1 = HE[nextH[he1]];
            int prev_e1 = HE[prevH[he1]];
            if (Emap[next_e1] != -1)
                c1 = int(!bool(Emap[next_e1]));
            else if (Emap[prev_e1] != -1)
                c1 = int(!bool(Emap[prev_e1]));

            // set color
            if (c0 >= 0 && c1 >= 0) {
                assert(c0 == c1); // make sure the colors agree
                return c0;
            }
            if (c0 >= 0) {
                return c0;
            }
            if (c1 >= 0) {
                return c1;
            }
        }
        return -1;
    };


    // --- at the end also go through all (including non-quad)faces and color any uncolored edges (do it twice, because the first time some edges might not be 'colorable', but then they are the second time)
    for (int attempt = 0; attempt < 2; ++attempt) {
        for (int fi = 0; fi < F.rows(); ++fi) {
            // if (D[fi] != 4) {
            for (int m = 0; m < D[fi]; ++m) {
                int ei = FE(fi, m);

                if (!is_boundary_edge(ei, EH)) { // if the edge is not on the boundary

                    // --- search using nextH
                    if (Emap[ei] == -1) { // if this edge it has not been colored
                        Emap[ei] = color_from_nextHs(ei, EH, twinH, nextH, HE, Emap);
                    }
                    else continue;

                    // --- if it was not found with the nextH, try using prevH instead
                    if (Emap[ei] == -1) {
                        Emap[ei] = color_from_prevHs(ei, EH, twinH, prevH, HE, Emap);
                    }
                    else continue;

                    // --- if it was still not found, try using the other edges of the neighboring triangles
                    if (Emap[ei] == -1) {
                        Emap[ei] = color_from_neighboring_edges(ei, EH, prevH, nextH, HE, Emap);
                    }
                    else continue;
                }
            }
            // }
        }
    }



    // --- make sure all boundary edges are uncolored
    for (int fi = 0; fi < F.rows(); ++fi)
        for (int m = 0; m < D[fi]; ++m) {
            int ei = FE(fi, m);
            if (is_boundary_edge(ei, EH))
                Emap[ei] = -1; // make sure it is uncolored
        }

    // --- check that all edges have been colored
    for (int ei = 0; ei < nE; ++ei)
        if (Emap[ei] == -1) {
            if (!is_boundary_edge(ei, EH)) {
                cerr << endl << "Attention! Edge ei = " << ei << " that is not on the boundary has NOT been colored (v0 = " << EV(ei, 0) << ". v1 = " << EV(ei, 1) << "). This shouldn't happen!" << endl << endl;
                //throw;
            }
        }

    // --- check if coloring is aligned with dominant direction, otherwise flip all entries of Emap
    if (!remember_existing_emap_values || e_remember_v1 == -1)
    {
        MatrixXd uv;
        uv.setZero(2, 2);

        int ei = 0; // pick an arbitrary edge that has color=0 
        while (ei < Emap.rows()) {
            if (Emap[ei] == 0) {
                int v1 = EV(ei, 0);
                int v2 = EV(ei, 1);
                uv.row(0) = UVcoords.row(v1);
                uv.row(1) = UVcoords.row(v2);
                if (abs(uv(0, 0) - uv(1, 0)) < 1.1 && abs(uv(0, 1) - uv(1, 1)) < 1.1) // make sure we are not on a cut where the UV jumps
                    break;
            }
            ++ei;
            if (ei == Emap.rows() - 1) { cerr << "While creating Emap, could not find any edge with color=0" << endl; throw; }
        }

        if (abs(uv(0, 0) - uv(1, 0)) > abs(uv(0, 1) - uv(1, 1)))  // if the subdominant direction changes more than the dominant
            MeshHelpers::flip_Emap(Emap);
    }
    else
    {
        // --- make sure Emap was maintained
        for (int ei = 0; ei < EV.rows(); ++ei) {
            int v1 = EV(ei, 0);
            int v2 = EV(ei, 1);
            if (min(v1, v2) == min(e_remember_v1, e_remember_v2) && max(v1, v2) == max(e_remember_v1, e_remember_v2)) // if we are on the same edge
                if (e_remember_Emap != Emap[ei])
                    MeshHelpers::flip_Emap(Emap);
        }
    }

}


void QuadMesh::add_extra_vertex_on_triangle_faces()
{
    // --- find pairs of triangles faces
    struct Pair {
        Pair(int f1, int f2, int ei) : ei(ei) {
            fis.push_back(f1);
            fis.push_back(f2);
        }
        vector<int> fis; // is filled with the two face indices of the pair
        int ei;
    };

    int nF = D.size();

    vector<Pair> pairs;
    vector<int> found_fis;
    for (int fi = 0; fi < nF; ++fi) {
        if (!Helpers::in_vector(found_fis, fi)) {
            if (D[fi] != 4) {
                for (int m = 0; m < D[fi]; ++m) {

                    int hei = FH(fi, m); // edge of fi
                    int ei = HE(hei);
                    assert(hei >= 0);


                    int twin = twinH[hei];
                    if (twin >= 0)
                    {
                        int nfi = HF(twin); // other neighboring face (contains the twin of hei)

                        if (!Helpers::in_vector(found_fis, nfi)) {
                            if (D[nfi] == D[fi] && D[fi] == 3) {
                                bool on_boundary = false; // make sure that fi and nfi are not on the boundary
                                if (is_boundary_face(fi, twinH, D, FH) || is_boundary_face(nfi, twinH, D, FH))
                                    on_boundary = true;

                                if (!on_boundary)
                                {
                                    pairs.emplace_back(fi, nfi, ei);
                                    found_fis.push_back(fi);
                                    found_fis.push_back(nfi);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "Found " << pairs.size() << " pairs of triangles: ";
    for (int i = 0; i < pairs.size(); ++i)  cout << " " << i << " : " << pairs[i].fis[0] << "," << pairs[i].fis[1] << endl;
    cout << endl;

    ////////////////////////////////////
    // --- add vertex on common edge

    for (const auto& pair : pairs) {
        int nV = V.rows(); // current number of vertices

        RowVector3d v = 0.5 * (V.row(EV(pair.ei, 0)) + V.row(EV(pair.ei, 1))); // new vertex coordinates
        V.conservativeResize(nV + 1, 3);
        V.row(nV) = v;

        // find the order on each face where the new vertex should be added
        int ei1 = EV(pair.ei, 0);
        int ei2 = EV(pair.ei, 1);

        vector<int> f_orders(2, -1); // is filled with 2 indices (1 for each face) : 0/1/2 if the new vertex should be added after the 0th/1st/2nd vertex respectively
        for (int k = 0; k < 2; ++k) { // for each of the two faces
            for (int i = 0; i < 3; ++i) {
                if ((F(pair.fis[k], i) == ei1 || F(pair.fis[k], i) == ei2) &&
                    (F(pair.fis[k], (i + 1) % 3) == ei1 || F(pair.fis[k], (i + 1) % 3) == ei2)) {
                    f_orders[k] = i;
                    break;
                }
            }
        }
        assert(f_orders[0] >= 0 && f_orders[1] >= 0); // make sure both were found

        for (int k = 0; k < 2; ++k) {
            int i = 3;
            while (i > f_orders[k])
            {
                if (i == f_orders[k] + 1)
                    F(pair.fis[k], i) = nV; // get new vi
                else
                    F(pair.fis[k], i) = F(pair.fis[k], i - 1); // get previous vi
                --i;
            }
        }

        D[pair.fis[0]] = 4;
        D[pair.fis[1]] = 4;
    }
}


void QuadMesh::get_UVcoords_on_quad_vertices(const directional::TriMesh& meshCut, const MatrixXd& triangle_mesh_UVcoords) {

    const MatrixXd& Vcut = meshCut.V;
    const MatrixXi& Fcut = meshCut.F;

    // --- get face centers of original (triangular) mesh
    const int nF = Fcut.rows();
    MatrixXd fcens;
    fcens.setZero(Fcut.rows(), 3);
    for (int fi = 0; fi < nF; ++fi)
        for (int m = 0; m < 3; ++m)
            fcens.row(fi) += Vcut.row(Fcut(fi, m)) / 3.0;


    // --- get UV coords for each Vquad from its closest original face
    const int nVquad = V.rows();

    UVcoords.setZero(nVquad, 2);
    VectorXi closest_fis(nVquad);
    for (int vi = 0; vi < nVquad; ++vi) {
        MatrixXd v = V.row(vi).replicate(nF, 1);
        VectorXd ds = (v - fcens).rowwise().squaredNorm();

        MatrixXf::Index cfi;
        ds.minCoeff(&cfi);

        //project Vquad.row(vi) on the found face (this step is actually not needed, because the proj is already on the face). Just doing it to check that the correct closest face was selected
        Vector3d N = meshCut.faceNormals.row(cfi); // normal of the plane of the face
        Vector3d proj = Geometry::project_pt_on_plane(N, Vector3d(fcens.row(cfi)), Vector3d(V.row(vi)));
        // cout << " RowVector3d(proj) - Vquad.row(vi)).squaredNorm() = " << (RowVector3d(proj) - Vquad.row(vi)).norm() << endl; 
        //assert((RowVector3d(proj) - Vquad.row(vi)).squaredNorm() < 1e-4);

        // get barycentric coordinates of Vquad on the closest original face
        const Vector3d B = MeshHelpers::barycentric_coordinates(Vcut.row(Fcut(cfi, 0)), Vcut.row(Fcut(cfi, 1)), Vcut.row(Fcut(cfi, 2)), proj);
        assert(0 <= B[0] <= 1);
        assert(0 <= B[1] <= 1);
        assert(0 <= B[2] <= 1);

        RowVector3d test(0, 0, 0);
        for (int m = 0; m < 3; ++m) {
            UVcoords.row(vi) += B[m] * triangle_mesh_UVcoords.row(Fcut(cfi, m));
            test += B[m] * Vcut.row(Fcut(cfi, m)); // check that B * V.row(v_original) = proj
        }
        //assert((test - RowVector3d(proj)).squaredNorm() < 1e-6);
        // cout << "vi = " << vi << " , uvs = " << np3dp->UVcoordsQuad.row(vi) << endl;
    }
}


MatrixXd QuadMesh::per_face_normals(const MatrixXd& V, const MatrixXi& T, const MatrixXi& F, const VectorXi& TF)
{
    MatrixXd N_triangles, N;
    VectorXd A_triangles;
    igl::per_face_normals(V, T, N_triangles);
    igl::doublearea(V, T, A_triangles);
    A_triangles *= 0.5;
    N.setZero(F.rows(), 3);
    for (int ti = 0; ti < T.rows(); ++ti) {
        int fi = TF[ti];
        N.row(fi) += N_triangles.row(ti) * A_triangles[ti];
    }
    N.rowwise().normalize();
    return N;
}


vector<vector<int>> QuadMesh::adjacency_list() const
{
    int nV = V.rows();
    vector<vector<int>> VV(nV);
    for (int ei = 0; ei < EV.rows(); ++ei) {
        int v0 = EV(ei, 0);
        int v1 = EV(ei, 1);
        VV[v0].push_back(v1);
        VV[v1].push_back(v0);
    }
	//Helpers::write_VV_to_txt(DATA_PATH + output_folder + "VV.txt", VV);
    return VV;
}



vector<vector<int>> QuadMesh::adjacency_list_with_Emap(int dir) const
{
	int nV = V.rows();
    vector<vector<int>> VV(nV);
    for (int ei = 0; ei < EV.rows(); ++ei) {
        if (Emap[ei] == dir)
        {
            int v0 = EV(ei, 0);
            int v1 = EV(ei, 1);
            VV[v0].push_back(v1);
            VV[v1].push_back(v0);
        }
    }
    return VV;
}


bool QuadMesh::shortest_path(int source_vertex, const set<int>& target_vertices, const vector<vector<int>>& adjacency_list, vector<int>& path) const
{
    VectorXd min_distance; VectorXi previous;
    igl::dijkstra(source_vertex, target_vertices, adjacency_list, min_distance, previous);

    // closest target
    int closest_target = -1; int min_dist = 1e4;
    for (int vi : target_vertices) {
        // cout << "vi : " << vi << " min dist = " << min_distance[vi] << endl;
        if (min_distance[vi] < min_dist) {
            closest_target = vi;
            min_dist = min_distance[vi];
        }
    }

    // shortest path
    int vi = closest_target;
    while (vi != source_vertex)
    {
        if (vi == -1) return false;
        path.push_back(vi);
        vi = previous[vi];
    }
    path.push_back(source_vertex);
    //cout << "path : " << endl;
    //Helpers::print_vector(path);

    return true;
}


int QuadMesh::count_Emap_values_of_face_edges(int fi, const VectorXi& D, const MatrixXi& FE, const VectorXi& Emap, int emap_value_to_count)
{
    int count = 0;
    for (int k = 0; k < D[fi]; ++k) {
        int ei = FE(fi, k);
        if (Emap[ei] == emap_value_to_count)
            ++count;
    }
    return count;
}



////////////////////////
/// --- Ribs on edges
////////////////////////

namespace RibsHelpers
{
    int get_next_ei_of_row(int ei, const MatrixXi& EH, const VectorXi& nextH, const VectorXi& prevH, const VectorXi& twinH, const VectorXi& HE)
    {
        int he = EH(ei, 0); if (he < 0) goto findnext;
        he = nextH[he];             if (he < 0) goto findnext;
        he = twinH[he];             if (he < 0) goto findnext;
        he = nextH[he];             if (he < 0) goto findnext;
        return HE[he];

    findnext:;
        he = EH(ei, 1); if (he < 0) return -1;
        he = prevH[he];         if (he < 0) return -1;
        he = twinH[he];         if (he < 0) return -1;
        he = prevH[he];         if (he < 0) return -1;
        return HE[he];
    }

    int get_prev_ei_of_row(int ei, const MatrixXi& EH, const VectorXi& nextH, const VectorXi& prevH, const VectorXi& twinH, const VectorXi& HE)
    {
        int he = EH(ei, 0); if (he < 0) goto findnext;
        he = prevH[he];             if (he < 0) goto findnext;
        he = twinH[he];             if (he < 0) goto findnext;
        he = prevH[he];             if (he < 0) goto findnext;
        return HE[he];

    findnext:;
        he = EH(ei, 1); if (he < 0) return -1;
        he = nextH[he];         if (he < 0) return -1;
        he = twinH[he];         if (he < 0) return -1;
        he = nextH[he];         if (he < 0) return -1;
        return HE[he];
    }

    int get_bridge_to_next_rib_ei(int dist, const MatrixXi& EH, int ei, const VectorXi& nextH, const VectorXi& twinH, const VectorXi& HE, const VectorXi& D, const VectorXi& HF, const VectorXi& ERibsMap)
    {
        int he = EH(ei, 0); //if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
        for (int m = 0; m < dist; ++m) {
            if (he < 0) return -1;
            he = nextH[he];
        	if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;  // if on boundary, or if crossing another rib-edge, or passing by a on a non-quad face
            he = nextH[he];
        	if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
            he = twinH[he];
        	if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
        }
        return HE[he];
    }

    int get_bridge_to_prev_rib_ei(int dist, const MatrixXi& EH, int ei, const VectorXi& prevH, const VectorXi& twinH, const VectorXi& HE, const VectorXi& D, const VectorXi& HF, const VectorXi& ERibsMap)
    {
        int he = EH(ei, 0); //if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
        for (int m = 0; m < dist; ++m) {
            if (he < 0) return -1;
            he = prevH[he];
            if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;   // if on boundary, or if crossing another rib-edge, or passing by a on a non-quad face
            he = prevH[he];
            if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
            he = twinH[he];
            if (he < 0) return -1; if (ERibsMap[HE[he]]) return -1;
        }
        return HE[he];
    }

    int min_dist_to_existing_rib(int ei, const VectorXi& ERibs, const MatrixXi& EH, const VectorXi& HF, const VectorXi& HE, const VectorXi& HV, const VectorXi& nextH, const VectorXi& twinH, const VectorXi& D) {
        int dist = 1e8;
        for (int m = 0; m < 2; ++m) {

            int current_step = 1;
            int he = EH(ei, 0);

            while (current_step < HF.size()) {
                if (he < 0 || D[HF[he]] != 4) { // conservatively also stop at boundaries or non-quad faces
                    // dist = min(dist, current_step);
                    break;
                }
                he = nextH[nextH[he]];
                if (ERibs[HE[he]] == 1) {
                    dist = min(dist, current_step);
                    break;
                }
                he = twinH[he]; if (he < 0) break;
                ++current_step;
            }
        }
        return dist;
    }

    vector<int> get_edge_loop(int start_ei, const VectorXi& Emap, const MatrixXi& EH, const VectorXi& nextH, const VectorXi& prevH, const VectorXi& twinH, const VectorXi& HE)
    {
        vector<int> eis;
        std::deque<int> Q = { start_ei };

        int emap_of_first = Emap[start_ei];

        while (!Q.empty()) {
            int ei = Q.front();
            Q.pop_front();

            const int next_ei = get_next_ei_of_row(ei, EH, nextH, prevH, twinH, HE);
            const int prev_ei = get_prev_ei_of_row(ei, EH, nextH, prevH, twinH, HE);

            if (next_ei >= 0) {
                if (Emap[next_ei] == emap_of_first) {
	                if (!Helpers::in_vector(eis, next_ei)) {
	                    eis.push_back(next_ei);
	                    Q.push_back(next_ei);
					}
                }

            }
            if (prev_ei >= 0) {
	            if (Emap[prev_ei] == emap_of_first) {
		            if (!Helpers::in_vector(eis, prev_ei)) {
		            	eis.push_back(prev_ei);
		            	Q.push_back(prev_ei);
		            }
	            }
            }
        }
        return eis;
    }
}


void QuadMesh::create_quad_ERibs()
{
    assert(FE.rows() == F.rows());
    int nE = EF.rows();

    ERibMap.setConstant(nE, 0); //
    // return;

    // --- select random sub-dominant edge and mark it as a rib
    int fi = ribs_f0;
    int start_ei = -1;
    while (start_ei < 0) {
	    for (int m=0; m<D[fi]; ++m) {
            const int ei = FE(fi, m);
            if (Emap[ei] == 1) {
	            start_ei = ei;
                break;
            }
	    }
        ++fi;
    }


	// --- get entire edge loop
    std::deque<int> Q = {start_ei};

    while(!Q.empty()) {
        int ei = Q.front();
        Q.pop_front();
        assert(Emap[ei] == 1);

        if (RibsHelpers::min_dist_to_existing_rib(ei, ERibMap, EH, HF, HE, HV, nextH, twinH, D) >= ribs_dist) {
            if (!ERibMap[ei]) {

                vector<int> ei_row = RibsHelpers::get_edge_loop(ei, Emap, EH, nextH, prevH, twinH, HE);

            	/*for (int ii : ei_row)
                	cout << EV(ii, 0) << endl << EV(ii, 1) << endl;
                int v0 = EV(ei, 0);
                int v1 = EV(ei, 1);
                cout << endl;*/

                bool too_close_to_existing_rib = false;
                for (int i : ei_row) {
                    if (RibsHelpers::min_dist_to_existing_rib(i, ERibMap, EH, HF, HE, HV, nextH, twinH, D) < ribs_dist) {
                        too_close_to_existing_rib = true;
                    }
                }

                if (!too_close_to_existing_rib) { // then set ei + ei_row to a rib 
                    ERibMap[ei] = 1;
                    for (int i : ei_row) {
                        if(Emap[i] != 1)
                        {
	                        cout << "Attention, the topology / Emap of the quad mesh is not correct. Cannot find ribs. Edge: ei = " << ei << ", vis : " << EV(ei, 0) << " , " << EV(ei, 1) << endl;
                        }// throw; }
                        ERibMap[i] = 1;

                        int other_ei = RibsHelpers::get_bridge_to_next_rib_ei(ribs_dist, EH, i, nextH, twinH, HE, D, HF, ERibMap);
                        if (other_ei >= 0)
                            if (!ERibMap[other_ei] && !Helpers::in_deque(Q, other_ei) && Emap[other_ei] == 1) // if it is not a rib, and it is not in Q, and it is a sub-dominant edge
                                Q.push_back(other_ei);

                        other_ei = RibsHelpers::get_bridge_to_prev_rib_ei(ribs_dist, EH, i, prevH, twinH, HE, D, HF, ERibMap);
                        if (other_ei >= 0)
                            if (!ERibMap[other_ei] && !Helpers::in_deque(Q, other_ei) && Emap[other_ei] == 1) // if it is not a rib, and it is not in Q, and it is a sub-dominant edge
                                Q.push_back(other_ei);
                    }
                }
            }
        }
    }

    /*for (int ei = 0; ei< ERibs.size(); ++ei)
        if (ERibs[ei])
            cout << EV(ei, 0) << endl << EV(ei, 1) << endl;
    cout << endl;*/
}


MatrixXd QuadMesh::get_ribs_polylines() const
{
    int n = ERibMap.sum();
    MatrixXd lines(n, 6);
    int i = 0;
    for (int ei = 0; ei < EV.rows(); ++ei) {
        if (ERibMap[ei]) {
            lines.block(i, 0, 1, 3) = V.row(EV(ei, 0));
            lines.block(i, 3, 1, 3) = V.row(EV(ei, 1));
            ++i;
        }
    }
    if (i != n) throw invalid_argument("Not all ribs edges were found"); // make sure we found all the ribs
    return lines;
}


///////////////////////////////////////////////
/// --- Functions related to partitioning
///////////////////////////////////////////////

void QuadMesh::expand_cutting_from_he(int he, VectorXi& C, vector<int>& passed_vis, int nE, const VectorXi& HE, const VectorXi& HV, const VectorXi& nextH, const VectorXi& twinH, const VectorXi& vertexDegree, const VectorXi& boundaryVertices, bool stop_at_cut_intersections, bool stop_at_irregular_vertices)
{
    if (C.rows() == 0) C.setZero(nE); // if C is not set, initialize it to 0
    else if (C.rows() != nE) throw invalid_argument("C.rows() != nE. Cannot expand cut from he");

    passed_vis.clear();
    do {
        //cout << HV[he] << endl;
        if (!Helpers::in_vector(passed_vis, HV[he])) passed_vis.push_back(HV[he]);
        if (he == -1) break;
        if (C[HE[he]] == 1) break; // if this edge has already been activated
        C[HE[he]] = 1; // activate cut of edge
        he = nextH[he];
        if (he == -1) break;
        if (stop_at_cut_intersections && C[HE[he]] == 1) break; // if crossing an edge that was already labeled as a cut, stop
        he = twinH[he];
        if (he == -1) break;
        he = nextH[he];
        if (he == -1) break;
        if (vertexDegree[HV[he]] != 4) break;
        if (twinH[he]!=-1) if (!Helpers::in_vector(passed_vis, HV[twinH[he]])) passed_vis.push_back(HV[twinH[he]]);
    } while ((!stop_at_irregular_vertices || vertexDegree[HV[he]] == 4) && boundaryVertices[HV[he]] != 1); // stop when you encounter another singularity, or a boundary
}

int QuadMesh::find_distance_to_boundary(int he) const
{
    vector<int> passed_eis;
    int count = 0;
    do {
        if (Helpers::in_vector(passed_eis, HE[he])) break;
        passed_eis.push_back(HE[he]);
        if (he == -1) break;
        he = nextH[he];
        if (he == -1) break;
        he = twinH[he];
        if (he == -1) break;
        he = nextH[he];
        if (he == -1) break;
        // cout <<  HV[he] << endl;
        ++count;
    } while (true); // stop when you encounter another singularity, or a boundary
    return count;
}


///////////////////////////////////////////////
/// --- Cleanup
///////////////////////////////////////////////

vector<int> QuadMesh::remove_boundary_isolated_triangles()
{
    cout << "Removing isolated boundary triangles... ";
    // no need for updated mesh data

    if (D.rows() != F.rows()) throw invalid_argument("Attention! D.rows() != F.rows()");
    
    MatrixXi newEV, newFE, newEF, newEFi; MatrixXd newFEs; VectorXi newinnerEdges;
    hedra::polygonal_edge_topology(D, F, newEV, newFE, newEF, newEFi, newFEs, newinnerEdges);

    //vector<int> vis_to_remove;
    vector<int> removed_fis;
    for (int fi = 0; fi < F.rows(); ++fi)
    {
        if (D[fi] == 3)
        {
            vector<int> eis;
            for (int k = 0; k < D[fi]; ++k)
            {
                int ei = newFE(fi, k);
                if (newEF(ei, 0) == -1 || newEF(ei, 1) == -1)
                    eis.push_back(ei);
            }
            if (eis.size() == 2)
            {
                //int common_vi;
                //if (newEV(eis[0], 0) == newEV(eis[1], 0) || newEV(eis[0], 0) == newEV(eis[1], 1))
                //    common_vi = newEV(eis[0], 0);
                //else
                //    common_vi = newEV(eis[0], 1);
                //vis_to_remove.push_back(common_vi);

                removed_fis.push_back(fi);
            }
        }
    }
    std::reverse(removed_fis.begin(), removed_fis.end());
    //std::sort(vis_to_remove.begin(), vis_to_remove.end(), [](const int& a, const int& b)->bool {return a > b; });

    for (int fi : removed_fis) {
        Helpers::removeRow(F, fi);
        Helpers::removeRow(D, fi);
    }
    //for (int vi : vis_to_remove)
    //{
    //    Helpers::removeRow(V, vi);
    //    Helpers::removeRow(UV, vi);
    //    for (int fi = 0; fi < F.rows(); ++fi)
    //        for (int k = 0; k < D[fi]; ++k)
    //            if (F(fi, k) >= vi)
    //                F(fi, k) -= 1;
    //}

    cout << "Completed. Removed " << removed_fis.size() << " faces." << endl;

    //vector<int> removed_vis = QuadMesh::cleanup_unreferenced_vertices(V, UV, FixedVis, F, D);
    return removed_fis;
}


vector<int> QuadMesh::remove_faces_without_dominant_edges()
{
    // Attention! This function requires an updated Emap!

    cout << "Removing faces without dominant edges... ";

    MatrixXi newEV, newFE, newEF, newEFi; MatrixXd newFEs; VectorXi newinnerEdges;
    hedra::polygonal_edge_topology(D, F, newEV, newFE, newEF, newEFi, newFEs, newinnerEdges);

    if (Emap.rows() != newEV.rows()) throw invalid_argument("Emap.rows()!= newEV.rows()");

    vector<int> removed_fis;
    for (int fi = 0; fi < F.rows(); ++fi)
    {
        int count_dominant = 0;
        for (int k = 0; k < D[fi]; ++k)
        {
            int ei = newFE(fi, k);
            if (Emap[ei] == 0)
                ++count_dominant;
        }
        if (count_dominant == 0)
            removed_fis.push_back(fi);
    }
    std::reverse(removed_fis.begin(), removed_fis.end());

    for (int fi : removed_fis) {
        Helpers::removeRow(F, fi);
        Helpers::removeRow(D, fi);
    }

    cout << "Completed. Removed " << removed_fis.size() << " faces." << endl;
    //vector<int> removed_vis = QuadMesh::cleanup_unreferenced_vertices(V, UV, FixedVis, F, D);
    return removed_fis;
}


vector<int> QuadMesh::cleanup_unreferenced_vertices()
{
    cout << "Cleaning up unreferenced vertices... ";

    if (V.rows()!= UVcoords.rows()) throw invalid_argument("V.rows()!= UV.rows()  " + to_string(V.rows()) + " != " + to_string(UVcoords.rows()));
    // if (V.rows()!= FixedVis.rows()) throw invalid_argument("V.rows()!= FixedVis.rows()  " + to_string(V.rows()) + " != " + to_string(FixedVis.rows()));

	VectorXi vis_encounters; vis_encounters.setZero(V.rows());

    for (int fi = 0; fi < F.rows(); ++fi) 
        for (int k = 0; k < D[fi]; ++k)
            ++vis_encounters[F(fi, k)];

    vector<int> unreferenced_vis;
    for (int vi=V.rows()-1; vi >=0; --vi) {
        if (vis_encounters[vi] == 0) {
            unreferenced_vis.push_back(vi);
            Helpers::removeRow(V, vi);
            Helpers::removeRow(UVcoords, vi);
            // Helpers::removeRow(FixedVis, vi);
            for (int fi = 0; fi < F.rows(); ++fi)
                for (int k = 0; k < D[fi]; ++k)
                    if (F(fi, k) >= vi)
                        F(fi, k) -= 1;
        }
    }

    cout << "Completed. Removed " << unreferenced_vis.size() << " unreferenced vis " << endl;

    return unreferenced_vis;
}


vector<int> QuadMesh::simplify_boundary_polyline()
{
    cout << "Simplifying boundary polylines... ";

    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
    hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);

    auto edge_is_on_boundary = [](int ei, const MatrixXi& EF)->bool { return EF(ei, 0) == -1 || EF(ei, 1) == -1; };

    vector<pair<int, int>> vis_fis_pairs_to_remove;
    for (int fi = 0; fi < F.rows(); ++fi)
    {
        if (D[fi] >= 4)
        {
            vector<int> boundary_eis;
            for (int k = 0; k < D[fi]; ++k)
            {
                int ei = FE(fi, k);
                if (edge_is_on_boundary(ei, EF)) {
	                boundary_eis.push_back(ei);
                }
            }

            if (boundary_eis.size() == 2)
            {
                const int ei1 = boundary_eis[0];
                const int ei2 = boundary_eis[1];

                // check if two edges have a common vertex
                int vi_common = -1;
                if (EV(ei1, 0) == EV(ei2, 0) || EV(ei1, 0) == EV(ei2, 1))
                    vi_common = EV(ei1, 0);
                else if (EV(ei1, 1) == EV(ei2, 0) || EV(ei1, 1) == EV(ei2, 1))
                    vi_common = EV(ei1, 1);
                if (vi_common == -1)
                    continue; // no vertex to collapse in this face

                // check angle between the two edges
                RowVector3d ei1_vec = (V.row(EV(ei1, 0)) - V.row(EV(ei1, 1))).normalized();
                RowVector3d ei2_vec = (V.row(EV(ei2, 0)) - V.row(EV(ei2, 1))).normalized();
                double dot = ei1_vec.dot(ei2_vec);
                double l1 = (V.row(EV(ei1, 0)) - V.row(EV(ei1, 1))).norm();
                double l2 = (V.row(EV(ei2, 0)) - V.row(EV(ei2, 1))).norm();

                if (abs(dot) < 0.75 && min(l1, l2)/max(l1,l2) > 0.05 && l1 > 1e-3 && l2 > 1e-3) 
                    continue; // too large angle could mean that it is a corner

                // store vi for removal
                vis_fis_pairs_to_remove.push_back(make_pair(vi_common, fi));
            }
        }
    }

    std::sort(vis_fis_pairs_to_remove.begin(), vis_fis_pairs_to_remove.end(), [](pair<int, int> fv1, pair<int, int> fv2) {return fv1.first > fv2.first; });

    vector<int> removed_vis;
    for (pair<int, int> fvi : vis_fis_pairs_to_remove)
    {
        int vi = fvi.first;
        int fi = fvi.second;
        int k = 0;
        for (k = 0; k < D[fi]; ++k)
            if (F(fi, k) == vi)
                break;

        // --- update vertex-based data
        Helpers::removeRow(V, vi);
        Helpers::removeRow(UVcoords, vi);
        // Helpers::removeRow(FixedVis, vi);
        removed_vis.push_back(vi);

        // --- update face-based data
        for (int j = k; j < D[fi]; ++j) { // first update current face
            if (j + 1 < F.cols())
                F(fi, j) = F(fi, j + 1);
            else
                F(fi, j) = -1;
        }

        for (int fj = 0; fj < F.rows(); ++fj) {// then update all faces
            for (int j = 0; j < D[fj]; ++j) {
                if (F(fj, j) == vi) throw invalid_argument("Attention! Cannot have this vertex on other faces. Fj = " + to_string(fj));
                if (F(fj, j) > vi)
                    F(fj, j) -= 1;
            }
        }

        D[fi] -= 1;
    }

    cout << "Completed. Removed " << vis_fis_pairs_to_remove.size() << " boundary vertices." << endl;

    return removed_vis;
}


vector<int> QuadMesh::remove_tiny_faces()
{
    double threshold = 1e-6;

    cout << "Removing tiny faces with area less than " << threshold << "... ";

    auto triangle_area = [](const RowVector3d& v1, const RowVector3d& v2, const RowVector3d& v3)->double
    {
        // get lengths of edges
        double a = (v1 - v2).norm();
        double b = (v2 - v3).norm();
        double c = (v3 - v1).norm();

        if (a < 0 || b < 0 || c < 0 || (a + b <= c) || a + c <= b || b + c <= a) {
            return 0;
        }
        float s = (a + b + c) / 2;
        return sqrt(s * (s - a) * (s - b) * (s - c));
    };

    vector<int> removed_fis;
    for (int fi = F.rows()-1; fi >=0; --fi)
    {
        double area = 0;

        if (D[fi] == 3) // find area of triangle
        {
            area = triangle_area(V.row(F(fi, 0)), V.row(F(fi, 1)), V.row(F(fi, 2)));
        }

        else if (D[fi] >3) // find area of polygon by summing the areas of the triangles
        {
            for (int j=0; j< D[fi]-2; ++j)
            {
                int vi1 = F(fi, 0);
                int vi2 = F(fi, j+1);
                int vi3 = F(fi, j+2);
	            area += triangle_area(V.row(vi1), V.row(vi2), V.row(vi3));
                if (area >= threshold)
                    break; // no need to calculate further
            }
        }

        if (area < threshold)
        {
	        Helpers::removeRow(F, fi);
            Helpers::removeRow(D, fi);
            removed_fis.push_back(fi);
        }
    }
    cout << "Completed! Removed " << removed_fis.size() << " faces." << endl;

    return removed_fis;
}


void QuadMesh::cleanup()
{
    cout << endl << "Cleaning up mesh..." << endl;

    //// --- remember edge-based manually selected partitions data
    //vector<vector<pair<int, int>>> remember_eis_vertices(2); // one for the manually selected eis of each side 
    //for (int side = 0; side < 2; ++side) {
    //    for (int ei : np3dp->manually_selected_partition_eis[side]) {
    //        int v0 = EV(ei, 0);
    //        int v1 = EV(ei, 1);
    //        remember_eis_vertices[side].push_back(make_pair(v0, v1));
    //    }
    //}

    if (V.rows() != UVcoords.rows()){ cerr << "V.rows()!= UVcoordsQuad.rows()  "  << V.rows() << " != " << UVcoords.rows() << endl; throw;}
    // if (V.rows() != np3dp->FixedVis.rows()) {cerr << "V.rows()!= FixedVis.rows()  " << V.rows() << " != "  << np3dp->FixedVis.rows() << endl; throw;}

    QuadMesh::update_quad_mesh_data();
    QuadMesh::create_quad_Emap(true);

    // --- Final cleanup steps
    const vector<int> removed_fis1 = QuadMesh::remove_faces_without_dominant_edges();
    const vector<int> removed_fis2 = QuadMesh::remove_tiny_faces();
    const vector<int> removed_fis3 = QuadMesh::remove_boundary_isolated_triangles();

	const vector<int> removed_vis1 = QuadMesh::simplify_boundary_polyline();
    const vector<int> removed_vis2 = QuadMesh::cleanup_unreferenced_vertices();

    // --- finally update all mesh data (if there have been updates)
    if (!removed_fis1.empty() || !removed_fis2.empty() || !removed_fis3.empty() || !removed_vis1.empty() || !removed_vis2.empty())
    {
        QuadMesh::update_quad_mesh_data();
        QuadMesh::create_quad_Emap(true);
    }

    // --- also update selection.selected_vis
	for (int& vi: selection.selected_vis) {
		for (int other_vi : removed_vis1)
            if (other_vi < vi)
                --vi;
        for (int other_vi : removed_vis2)
            if (other_vi < vi)
                --vi;
	}

    // --- also update np3dp->vis_to_check_for_alignments
    for (pair<int, int>& vis : vis_to_check_for_alignments) {
        for (int other_vi : removed_vis1){
            if (other_vi < vis.first)
                --vis.first;
            if (other_vi < vis.second)
                --vis.second;
        }
        for (int other_vi : removed_vis2){
            if (other_vi < vis.first)
                --vis.first;
            if (other_vi < vis.second)
                --vis.second;
        }
    }


    //// --- restore edge-based manually selected partitions data
    //vector<vector<int>> removed_vis = { removed_fis_vis2.second, removed_fis_vis1.second, removed_vis1, removed_fis_vis3.second, removed_vis2};
    //restore_saved_edges(np3dp, remember_eis_vertices, removed_vis);

    cout << "Cleaning up completed." << endl << endl;

    if (V.rows() != UVcoords.rows()) {cerr << "V.rows()!= UVcoordsQuad.rows()  " << V.rows() << " != " << UVcoords.rows() << endl; throw;}
}


void QuadMesh::save_data() const
{
    string data_folder = DATA_PATH + output_folder;
    igl::writeOBJ(data_folder + "quad_mesh.obj", V, F); cout << "Saved : " << data_folder + "quad_mesh.obj" << endl;
    hedra::polygonal_write_OFF(data_folder + "quad_mesh.off", V, D, F); cout << "Saved to file : " << data_folder + "quad_mesh.off" << endl;
    Helpers::write_matrix_to_txt(V, data_folder + "V.txt");

    MatrixXi Fcopy = F; // replace 0s with -1s
    for (int fi=0; fi<F.rows(); ++fi){
	    for (int k = D[fi]; k<F.cols(); ++k){
            Fcopy(fi, k) = -1;
	    }
    }
    Helpers::write_matrix_to_txt(Fcopy, data_folder + "F.txt");
    Helpers::write_Emap_to_txt(data_folder + "Emap.txt", Emap, EV);
}


///////////////////////////////////
/// --- Collapse
///////////////////////////////////

void QuadMesh::collapse_strip(VectorXi& strip, int selected_strip_direction)
{
    cout << "Collapsing selected strip...";

    int other_direction = static_cast<int>(!static_cast<bool>(selected_strip_direction));

    update_quad_mesh_data();
    if (!has_updated_mesh_data()) {cerr << "QuadMeshEditing::collapse_strip() requires updated mesh data" << endl; throw;}

    remember_previous_data();

    // saving original quad mesh (for debugging)
    //MeshHelpers::write_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, DATA_PATH + np3dp->output_folder);

    // find vertices that are in the center of a strip self-folding
    VectorXi immovable_vis; immovable_vis.setZero(V.rows());

    //cout << "Immovable vis : " << endl; // find vis who are in the middle of a loop
    for (int ei = 0; ei < EV.rows(); ++ei) {
        if (Emap[ei] != other_direction) { // only consider edges in the strip direction, (in the other direction all edges have both faces in the strip)
            int f1 = EF(ei, 0);
            int f2 = EF(ei, 1);
            if (f1 >= 0 && f2 >= 0) {
                if (strip[f1] == 1 && strip[f2] == 1) {
                    int v1 = EV(ei, 0);
                    int v2 = EV(ei, 1);
                    if (vertexDegree[v1] > 2 || boundaryVertices[v1]) { // accept all regular + boundary vertices 
                        immovable_vis[v1] = 1;
                        //cout << v1 << endl;
                    }
                    if (vertexDegree[v2] > 2 || boundaryVertices[v2]) { // accept all regular + boundary vertices 
                        immovable_vis[v2] = 1;
                        //cout << v2 << endl;
                    }
                }
            }
        }
    }


    // find which vertices need to be collapsed
    MatrixXi e_merge_plans;
    e_merge_plans.setConstant(EV.rows(), 3, -1); // 1st col: vi, 2nd col: vj, 3rd col: 0/1/-1

    for (int ei = 0; ei < EV.rows(); ++ei) {
        int f1 = EF(ei, 0);
        int f2 = EF(ei, 1);

        // --- decide if edge should be collapsed
        bool to_collapse = true;
	    {
        	// only one of the two Emap directions should be collapsed
            if (Emap[ei] == selected_strip_direction) to_collapse = false;
		    
		    //for non-boundary edges 
        	if (f1 >= 0) if (strip[f1] == 0) to_collapse = false;
        	if (f2 >= 0) if (strip[f2] == 0) to_collapse = false;

        	// for boundary edges of strip
            if (to_collapse && (f1 == -1 || f2 == -1))
            {
                int f = f1 == -1 ? f2 : f1; // find face that is not a boundary
                if (D[f] == 4 && count_Emap_values_of_face_edges(f, D, FE, Emap, other_direction) == 2) 
                    to_collapse = false; // if D[fi] = 4: don't collapse. Otherwise (for d=5) collapse
                else if (D[f] == 3) 
                    to_collapse = false;
            }
	    }

        if (to_collapse) 
        {
            int v1 = EV(ei, 0);
            int v2 = EV(ei, 1);

            if (v1 == v2) continue; // this should never happen, not sure why it comes up...

            e_merge_plans(ei, 0) = v1;
            e_merge_plans(ei, 1) = v2;

            if ((immovable_vis[v1] && !immovable_vis[v2]) || (boundaryVertices[v1] && !boundaryVertices[v2])) { // --- merge in location of v1
                e_merge_plans(ei, 2) = 0; // merge v2 on v1
            }
            else if ((immovable_vis[v2] && !immovable_vis[v1]) || (boundaryVertices[v2] && !boundaryVertices[v1])) { // --- merge in location of v2
                e_merge_plans(ei, 2) = 1; // merge v1 on v2
            }
            else { // ---  merge at the middle
                e_merge_plans(ei, 2) = -1; // merge v1 + v2 in the middle
            }
        }
    }

    // collapse
    int nV = V.rows(); // original number of vertices
    VectorXi V_to_Vcollapsed(nV); for (int vi = 0; vi < nV; ++vi) V_to_Vcollapsed[vi] = vi;
    MatrixXi newF(F.rows() - strip.sum(), F.cols());
    VectorXi newD(newF.rows());

    // --- execute e plans by merging vertices
    for (int ei = 0; ei < e_merge_plans.rows(); ++ei)
    {
        if (e_merge_plans(ei, 0) != -1)
        {
            const int v1_original = e_merge_plans(ei, 0);
            const int v2_original = e_merge_plans(ei, 1);
            const int v1 = V_to_Vcollapsed[v1_original];
            const int v2 = V_to_Vcollapsed[v2_original];
            const int type = e_merge_plans(ei, 2);

            if (v1 == v2) continue; // skip this if two vertices are already identical (this should never happen...)

            // remember V and UV data of v1 and v2
            RowVector3d V1 = V.row(v1);
            RowVector3d V2 = V.row(v2);
            RowVector2d UV1 = UVcoords.row(v1);
            RowVector2d UV2 = UVcoords.row(v2);
            // bool fixed = np3dp->FixedVis[v1] == 1 && np3dp->FixedVis[v2] == 1;

            // remove v1 and v2
            Helpers::removeRow(V, max(v1, v2));
            Helpers::removeRow(V, min(v1, v2));
            Helpers::removeRow(UVcoords, max(v1, v2));
            Helpers::removeRow(UVcoords, min(v1, v2));
            // Helpers::removeRow(np3dp->FixedVis, max(v1, v2));
            // Helpers::removeRow(np3dp->FixedVis, min(v1, v2));

            // remember which old V_to_Vcollapsed were pointing to v1 or v2 (these should also be updated to point to nvi)
            vector<int> vis_pointing_to_collapsed;
            for (int vi = 0; vi < nV; ++vi)
                if (vi != v1_original && vi != v2_original)
                    if (V_to_Vcollapsed[vi] == v1 || V_to_Vcollapsed[vi] == v2)
                        vis_pointing_to_collapsed.push_back(vi);

            // update V_to_Vcollapsed (relevant for all vertices that are not v1_original, v2_original or  )
            for (int vi = 0; vi < nV; ++vi) {
                if (min(v1, v2) < V_to_Vcollapsed[vi] && V_to_Vcollapsed[vi] < max(v1, v2))
                    V_to_Vcollapsed[vi] -= 1;
                else if (max(v1, v2) < V_to_Vcollapsed[vi])
                    V_to_Vcollapsed[vi] -= 2;
            }

            // create new vertex
            int nvi = V.rows();
            V.conservativeResize(nvi + 1, 3);
            UVcoords.conservativeResize(nvi + 1, 2);
            // np3dp->FixedVis.conservativeResize(nvi + 1); np3dp->FixedVis[nvi] = fixed ? 1 : 0;

            // update V_to_Vcollapsed
            V_to_Vcollapsed[v1_original] = nvi;
            V_to_Vcollapsed[v2_original] = nvi;
            for (int vi : vis_pointing_to_collapsed) 
                V_to_Vcollapsed[vi] = nvi;

            // fill in new entry in V and UV
            if (type == 0) {
                V.row(nvi) = V1;
                UVcoords.row(nvi) = UV1;
            }
            else if (type == 1) {
                V.row(nvi) = V2;
                UVcoords.row(nvi) = UV2;
            }
            else {
                V.row(nvi) = 0.5 * (V1 + V2);
                UVcoords.row(nvi) = 0.5 * (UV1 + UV2);
            }

            //Helpers::write_matrix_to_txt(np3dp->Vquad, DATA_PATH + np3dp->output_folder + "V.txt");
            //cout << endl;
        }
    }

    // collapse faces
    int fi_count = 0;
    for (int fi = 0; fi < F.rows(); ++fi) {
        if (strip[fi] == 0) {
            newF.row(fi_count) = F.row(fi);
            for (int k = 0; k < D[fi]; ++k) {
                int vi = newF(fi_count, k);
                newF(fi_count, k) = V_to_Vcollapsed[vi]; // update vertex indices
            }
            newD[fi_count] = D[fi];
            ++fi_count;
        }
    }
    if (fi_count != newF.rows()) {cerr << "Error while collapsing strip: fi_count = " << fi_count << " , newF.rows() = "  << newF.rows() << endl; throw;}

    strip.setZero(F.rows()); // reset strip

    // Update np3dp->selected_vis (if they are filled in)
    for (int& vi : selection.selected_vis)
        vi = V_to_Vcollapsed[vi];

    // Update np3dp->vis_to_check_for_alignments
    for (pair<int, int>& vis : vis_to_check_for_alignments){
        vis.first = V_to_Vcollapsed[vis.first];
        vis.second = V_to_Vcollapsed[vis.second];
    }
    cout << "Completed." << endl;

    // ------------------------------------------------------------------------------ Cleaning up --->
    // --- remove identical vertices from faces
    for (int fi = newF.rows() - 1; fi >= 0; --fi) {
        vector<int> passed_vis;
        for (int k = 0; k < newD[fi]; ++k) {
            int vi = newF(fi, k);
            if (Helpers::in_vector(passed_vis, vi)) { // if there's a duplicate vertex
                if (newD[fi] == 3) { // then remove face completely
                    Helpers::removeRow(newF, fi);
                    Helpers::removeRow(newD, fi);
                    goto ctn;
                }
                if (newD[fi] > 3) // then just remove that vertex
                {
                    do {
                        if (k < newF.cols() - 1)
                            newF(fi, k) = newF(fi, k + 1);
                        else
                            newF(fi, k) = -1;
                        ++k;
                    } while (k < newD[fi]);
                    newD[fi] -= 1;
                    goto ctn;
                }
            }
            passed_vis.push_back(vi);
        }
    ctn:;
    }

    // replace F, D with newF, newD
    F= newF;
    D = newD;

    if (V.rows() != UVcoords.rows()) {cerr << "Error while collapsing strip: V.rows()!= UV.rows()  "  << V.rows() << " != " + to_string(UVcoords.rows()) << endl; throw; }
    // if (V.rows() != np3dp->FixedVis.rows()) {cerr << "Error while collapsing strip: V.rows()!= FixedVis.rows()  " << V.rows() << " != " << np3dp->FixedVis.rows() << endl; throw;}

    cleanup();
    save_data();

    if (!has_updated_mesh_data()) throw invalid_argument("QuadMeshEditing::collapse_strip() finished without updated mesh data");
}


void QuadMesh::collapse_edge(int ei)
{
    // remember previous data (enables Undo operation)
    remember_previous_data();

    cout << "Collapsing edge ei : " << ei << " with vis : " << EV(ei, 0) << " , " << EV(ei, 1) << endl;
    int v1 = min(EV(ei, 0), EV(ei, 1));
    int v2 = max(EV(ei, 0), EV(ei, 1));
    int f1 = EF(ei, 0);
    int f2 = EF(ei, 1);

    for (int fi : {f1, f2}) { // if v1 and v2 are the first and last vertex of the face, then shift all vertices by one position, otherwise the vertex removal won't work
        if (F(fi, 0) == v1 && F(fi, D[fi] - 1) == v2 || F(fi, 0) == v2 && F(fi, D[fi] - 1) == v1) {
            int last_vi = F(fi, D[fi] - 1);
            for (int k = D[fi] - 1; k > 0; --k)
                F(fi, k) = F(fi, k - 1);
            F(fi, 0) = last_vi;
        }
    }

    // vertices
    V.row(v1) = 0.5 * (V.row(v1) + V.row(v2));
    Helpers::removeRow(V, v2);

    // UV
    UVcoords.row(v1) = 0.5 * (UVcoords.row(v1) + UVcoords.row(v2));
    Helpers::removeRow(UVcoords, v2);

    // FixedVis
    // Helpers::removeRow(np3dp->FixedVis, v2);

    // update np3dp->selected_vis
    for (int& vi : selection.selected_vis) {
        if (v1 < vi) --vi;
		if (v2 < vi) --vi;
    }

    // update np3dp->vis_to_check_for_alignments
    for (pair<int, int>& vis : vis_to_check_for_alignments){
        if (v1 < vis.first)  --vis.first;
        if (v1 < vis.second) --vis.second;
        if (v2 < vis.first)  --vis.first;
        if (v2 < vis.second) --vis.second;
    }

    // faces
    for (int fi = 0; fi < F.rows(); ++fi) {
        if (fi == f1 || fi == f2) {
            bool vi_shifting = false;
            for (int k = 0; k < D[fi]; ++k) {
                if (vi_shifting) {
                    if (k < F.cols() - 1)
                        F(fi, k) = F(fi, k + 1);
                    else
                        F(fi, k) = -1;
                }
                else if (F(fi, k) == v1 || F(fi, k) == v2) {
                    vi_shifting = true;
                    F(fi, k) = v1;
                }
            }
            assert(vi_shifting);
        }

        for (int k = 0; k < D[fi]; ++k)
            if (F(fi, k) == v2)
                F(fi, k) = v1;
            else if (F(fi, k) > v2)
                F(fi, k) -= 1;
    }

    // D
    D[f1] -= 1;
    D[f2] -= 1;
    
    // update mesh data
    update_quad_mesh_data();
    create_quad_Emap();
    save_data();
}


///////////////////////////////////
/// --- Smoothing 
///////////////////////////////////

namespace Smoothen_helpers
{
	void extend_boundary_edges(double mult, const VectorXi& prevH, const VectorXi& twinH, const VectorXi& nextH, const VectorXi& VH, const MatrixXi& EV,
        const VectorXi& HV, const VectorXi& HE, VectorXi& is_Vi_modified, VectorXd& V_ts, vector<vector<int>>& adjacent_non_boundary_edges, MatrixXd& V)
	{
        const int number_of_references = 2;


        is_Vi_modified.setZero(V.rows());
        V_ts.setZero(V.rows());
        for (int vi = 0; vi < V.rows(); ++vi)
        {
            if (adjacent_non_boundary_edges[vi].size() == 1)
            {
                int ei = adjacent_non_boundary_edges[vi][0];
                int other_vi = EV(ei, 0) == vi ? EV(ei, 1) : EV(ei, 0);

                // --- find direction
                RowVector3d dir = V.row(vi) - V.row(other_vi);
                double len = dir.norm();
                dir.normalize();

                // --- find reference edges
                vector<int> ei_refs;
                //int ei_ref = -1;

                int he = VH(other_vi, 0); // get halfedge going from other_vi to vi
                while (HV[twinH[he]] != vi && he != -1)
                    he = nextH[twinH[he]];

                int count = 0;
                while (he != -1 && count < number_of_references)
                {
                    he = prevH[he];
                    if (he != -1) he = twinH[he];
                    if (he != -1) he = prevH[he];
                    if (he != -1) ei_refs.push_back(HE[he]); //ei_ref = HE[he];
                    ++count;
                }

                if (ei_refs.empty()) continue;

                // --- get reference direction and length
                RowVector3d dir_ref(0, 0, 0);
                double len_ref = 0;

                for (int ei_ref : ei_refs)
                {
                    int vi_ref_0 = EV(ei_ref, 0);
                    int vi_ref_1 = EV(ei_ref, 1);
                    len_ref += (V.row(EV(ei_ref, 0)) - V.row(EV(ei_ref, 1))).norm();
                    RowVector3d current_dir = (V.row(vi_ref_0) - V.row(vi_ref_1)).normalized();
                    if (current_dir.dot(dir) < 0)
                        current_dir *= -1;
                    dir_ref += current_dir;
                }
                dir_ref.normalize();
                len_ref /= static_cast<double>(ei_refs.size());

                double t = len / len_ref;

                // --- move vertex to 'extended' position
                V.row(vi) = V.row(other_vi) + dir_ref * len_ref * mult;

                is_Vi_modified[vi] = 1;
                V_ts[vi] = t;
            }
        }
	}

    void cut_boundary_edges_to_original(double mult, const MatrixXi& EV, const VectorXi& is_Vi_modified, const VectorXd& V_ts, const vector<vector<int>>& adjacent_non_boundary_edges, MatrixXd& V)
	{
        for (int vi = 0; vi < V.rows(); ++vi)
        {
            if (is_Vi_modified[vi])
            {
                int ei = adjacent_non_boundary_edges[vi][0];
                int other_vi = EV(ei, 0) == vi ? EV(ei, 1) : EV(ei, 0);

                RowVector3d dir = V.row(vi) - V.row(other_vi);
                double len = dir.norm();
                dir.normalize();

                V.row(vi) = V.row(other_vi) + dir * len * V_ts[vi] / mult;
            }
        }
	}
}



void QuadMesh::smoothen(const directional::TriMesh& tri_mesh, bool with_smoothened_boundary, int iterations, bool with_cleanup)
{
    if (!has_updated_mesh_data()) throw invalid_argument("QuadMeshEditing::smoothen(Np3dp& np3dp, bool with_boundaries) requires updated mesh data");

    cout << "Smoothing with " << iterations << " iterations" << endl;

	//debug
	//MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);

    // remember previous data (enables Undo operation)
    remember_previous_data();

    const MatrixXd Voriginal = V; // create a copy

    int nV = V.rows();

    // --- find vertices on boundary and original boundary polyline
    MatrixXd boundaryPolyline; VectorXi boundary_vis;
    MeshHelpers::get_boundary_polyline_and_vertices(MatrixXd(), boundary_vis, V, T);
    MeshHelpers::get_boundary_polyline_and_vertices(boundaryPolyline, VectorXi(), tri_mesh.V, tri_mesh.F); // get boundary polyline from original triangle mesh

    // --- get avg edge length
    double avg_e_len = 0.0;
    for (int ei=0; ei<EV.rows(); ++ei){
	    int v1 = EV(ei, 0);
        int v2 = EV(ei, 1);
        avg_e_len += (V.row(v1) - V.row(v2)).norm();
    }
    avg_e_len /= static_cast<double>(EV.rows());


    // --- create mask coefficients
    VectorXd coefficients; coefficients.setOnes(nV); // set dummy values to coefficients
    if (with_smoothened_boundary) { // then (fake) boundary vertices are completely fixed, all other vertices move freely
        for (int vi=0; vi<V.rows(); ++vi)
            if (boundary_vis[vi])
                coefficients[vi] = 0;
    }
    else // then create smooth 'mask' coefficients where vertices close to the boundary move less and less.
    {
        if (boundaryPolyline.rows() > 0)
        {
            VectorXd bndry_dist(nV);
            for (int vi = 0; vi < nV; ++vi)
            {
                Vector3d pt = V.row(vi);
                bndry_dist(vi) = (MeshHelpers::projectToPolyline(pt, boundaryPolyline) - pt).norm();
            }
            double threshold = 3.5 * avg_e_len;
            for (int vi = 0; vi < nV; ++vi) {
                if (bndry_dist(vi) < threshold){
                    coefficients[vi] = Helpers::remap_unbound(bndry_dist(vi), 0.0, threshold, 0.0, 1.0); // coefficients go from 1.0 to 0.0 as we approach the boundary
                    coefficients[vi] *= coefficients[vi]; // squared
                }
                else{
                    coefficients[vi] = 1.0;
                }
                    
            }
        }
    }
	Helpers::write_vector_to_txt(coefficients, DATA_PATH + output_folder + string("coefficients.txt"));


    // --- create adjacency lists
    vector<vector<int>> Adjacency(nV);
    vector<vector<int>> Boundary_adjacency(nV);

    auto edge_is_on_boundary = [](int ei, const MatrixXi& EF)->bool { if (EF(ei, 0) == -1 || EF(ei, 1) == -1) return true; return false; };

    for (int ei = EV.rows()-1; ei >=0 ; --ei) {
	        int v1 = EV(ei, 0);
	        int v2 = EV(ei, 1);
        if (!edge_is_on_boundary(ei, EF)) { 
    		Adjacency[v1].push_back(v2);
    		Adjacency[v2].push_back(v1);
        }
        else {
    		Boundary_adjacency[v1].push_back(v2);
    		Boundary_adjacency[v2].push_back(v1);
        }
    }

    // --- debug --->
	//Helpers::write_VV_to_txt(DATA_PATH + np3dp->output_folder + "VV.txt", Adjacency);
	//Helpers::write_VV_to_txt(DATA_PATH + np3dp->output_folder + "VV_boundary.txt", Boundary_adjacency);
    // <---


    // --- construct uniform Laplacians (using adjacency lists) : nVxnV
    SparseMatrix<double> L;
    {
		vector<Triplet<double>> L_triplets;
		for (int vi = 0; vi < nV; ++vi)
		{
		    // off diagonal
		    for (int nv : Adjacency[vi])
		        L_triplets.emplace_back(vi, nv, 1.0 / static_cast<double>(Adjacency[vi].size()));
		    // on the diagonal
		    L_triplets.emplace_back(vi, vi, -1.0);
		}
		L.conservativeResize(nV, nV);
		L.setFromTriplets(L_triplets.begin(), L_triplets.end());
    }

	SparseMatrix<double> L_boundary;
    {
		vector<Triplet<double>> L_triplets;
		for (int vi = 0; vi < nV; ++vi)
		{
            if (!Boundary_adjacency[vi].empty())
            {
			    // off diagonal
			    for (int nv : Boundary_adjacency[vi])
			    {
			        L_triplets.emplace_back(vi, nv, 1.0 / static_cast<double>(Boundary_adjacency[vi].size()));
			    }
				// on the diagonal
				L_triplets.emplace_back(vi, vi, -1.0);
            }
		}
		L_boundary.conservativeResize(nV, nV);
		L_boundary.setFromTriplets(L_triplets.begin(), L_triplets.end());
    }

    
    // --- get per-vertex mass matrix (VxV)
    VectorXd A; // per face areas of triangles
    igl::doublearea(V, T, A);
    A *= 0.5;
    A = MeshHelpers::FtoV(V, T) * A; // convert to per-vertex areas
    SparseMatrix<double> M;
    Helpers::diag(A, M);

    SparseMatrix<double> Id(nV, nV); Id.setIdentity();

    // --- get delta step size find the smallest non-zero eigenvalue of the system L* X = s * M * X 
    double delta = 0.01; // default value in case the eigenvalue problem doesnt converge
    double boundary_delta = 0.5;
	{
        // make positive definite
    	L = L * -1.0; 
        L_boundary = L_boundary * -1.0;

        // --- find smallest positive eigenvalue using matlab
        const auto get_step_size = [](const SparseMatrix<double>& M, const SparseMatrix<double>& L, double default_value)->double {

            Engine* engine;

            MatrixXd EV; // Eigenvectors of the laplacian (w. mass matrix)
            MatrixXd ED; // Eigenvalues of the laplacian (w. mass matrix)

            // Launch MATLAB
            igl::matlab::mlinit(&engine);

            // Send Laplacian matrix to matlab
            igl::matlab::mlsetmatrix(&engine, "L", L);

            // Send mass matrix to matlab
            igl::matlab::mlsetmatrix(&engine, "M", M);

            // Extract the first 5 eigenvectors
            igl::matlab::mleval(&engine, "[EV,ED] = eigs(L,M,5,'smallestabs')");
            // Turn eigenvalue diagonal matrix into a vector
            igl::matlab::mleval(&engine, "ED=diag(ED)");

            // Retrieve the result
            igl::matlab::mlgetmatrix(&engine, "EV", EV);
            igl::matlab::mlgetmatrix(&engine, "ED", ED);

            for (int i = 0; i < ED.rows(); ++i)
            {
                if (ED(i, 0) > 1e-6)
                {
                    return 0.001 / ED(i, 0);
                }
            }
            return default_value; // just return default value
        };

        delta = get_step_size(M, L, delta);
        // boundary_delta = get_step_size(Id, L_boundary, boundary_delta);
        cout << "Step sizes: delta = " << delta << ", boundary_delta = " << boundary_delta << endl;

        // return to original sign
        L = L * -1.0; 
        L_boundary = L_boundary * -1.0;
	}
    

    // --- smoothening iterations
    cout << "--- Smoothing quad mesh with " << iterations << " iterations, step size = " << delta << ", and projecting it back on original surface." << endl;

	//// --- explicit steps
	//for (int iteration = 0; iteration < iterations; ++iteration) {
	//   V = (Id + delta * coefficients.asDiagonal() * L)* V;
	//}

    // --- implicit steps
	SparseLU<SparseMatrix<double>> solver(M - delta * coefficients.asDiagonal() * L); //Id - delta * coefficients.asDiagonal() * L
    if (solver.info() != Eigen::Success) { cerr << "Could not compute implicit solver. Exiting without smoothing." << endl; return; }

    SparseLU<SparseMatrix<double>> solver_boundary;
	if (with_smoothened_boundary){
        solver_boundary.compute(Id - boundary_delta * L_boundary); //Id - delta * coefficients.asDiagonal() * L
        if (solver_boundary.info() != Eigen::Success) { cerr << "Could not compute implicit solver. Exiting without smoothing." << endl; return; }
    }

	for (int iteration = 0; iteration < iterations; ++iteration) {
        cout << "Smoothening iteration " << iteration << endl;

        if (with_smoothened_boundary) // first do laplacian smoothing on the boundary only
        {
			for (int i=0; i<5; ++i)
			{
                // smoothen
                V = solver_boundary.solve(V).eval();

                // project to original boundary
				for (int vi = 0; vi < nV; ++vi) 
					if (boundary_vis[vi] == 1)
						V.row(vi) = MeshHelpers::projectToPolyline(V.row(vi), boundaryPolyline);
			}
            //MeshHelpers::save_mesh_data(V, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);
        }

	    V = solver.solve(M * V).eval();
		//MeshHelpers::save_mesh_data(V, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);

	    // --- project to original mesh surface + boundary
	    VectorXd sqrD; // #P list of smallest squared distances
	    VectorXi I; // #P list of primitive indices corresponding to smallest distances
	    MatrixXd C; // #P by 3 list of closest points
	    //igl::point_mesh_squared_distance(V, Voriginal, T, sqrD, I, C);
	    igl::point_mesh_squared_distance(V, tri_mesh.V, tri_mesh.F, sqrD, I, C);
	    V = C;
        //MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);


	    // --- project boundary to original boundary
        {
            for (int vi = 0; vi < nV; ++vi)
                if (boundary_vis[vi] == 1)
                    V.row(vi) = MeshHelpers::projectToPolyline(V.row(vi), boundaryPolyline);
        }


        //MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);
    }

    // ----------------------------- Cleanup ----------------------------------------------------------------------------->
    if (with_cleanup)
    {
        const vector<int> removed_vis1 = simplify_boundary_polyline();
        const vector<int> removed_fis = remove_tiny_faces();
        const vector<int> removed_vis2 = cleanup_unreferenced_vertices();

        if (!removed_vis1.empty() || !removed_fis.empty() || !removed_vis2.empty())
        {
            update_quad_mesh_data();
            create_quad_Emap();
        }

        // --- also update np3dp->selected_vis
        for (int& vi : selection.selected_vis) {
            for (int other_vi : removed_vis1)
                if (other_vi < vi)
                    --vi;
            for (int other_vi : removed_vis2)
                if (other_vi < vi)
                    --vi;
        }

        // --- also update np3dp->vis_to_check_for_alignments
        for (pair<int, int>& vis : vis_to_check_for_alignments) {
            for (int other_vi : removed_vis1) {
                if (other_vi < vis.first)
                    --vis.first;
                if (other_vi < vis.second)
                    --vis.second;
            }
            for (int other_vi : removed_vis2) {
                if (other_vi < vis.first)
                    --vis.first;
                if (other_vi < vis.second)
                    --vis.second;
            }
        }
    }
    save_data();

    if (!QuadMesh::has_updated_mesh_data()) {cerr << "Attention! QuadMeshEditing::smoothen(Np3dp& np3dp, bool with_boundaries) finished without updated mesh data" << endl; throw;}
}


void QuadMesh::smoothen_specific_quad_mesh(MatrixXd& V, const MatrixXi& F, const VectorXi& D, const MatrixXd& Triangle_V, const MatrixXi& Triangle_F, bool with_smoothened_boundary, int iterations)
{

    cout << "Smoothing with " << iterations << " iterations" << endl;

	//debug
	//MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);

    MatrixXi T; VectorXi TF;
    hedra::triangulate_mesh(D, F, T, TF);

    const MatrixXd Voriginal = V;

    int nV = V.rows();

    // --- get halfedge datastructure
    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
    MatrixXi EH, FH; VectorXi VH, HE, HF, HV, nextH, prevH, twinH;
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

    // --- find vertices on boundary and original boundary polyline
    MatrixXd boundaryPolyline; VectorXi boundary_vis;
    MeshHelpers::get_boundary_polyline_and_vertices(MatrixXd(), boundary_vis, V, T);
    MeshHelpers::get_boundary_polyline_and_vertices(boundaryPolyline, VectorXi(), Triangle_V, Triangle_F); // get boundary polyline from original triangle mesh

    // --- get avg edge length
    double avg_e_len = 0.0;
    for (int ei=0; ei<EV.rows(); ++ei){
	    int v1 = EV(ei, 0);
        int v2 = EV(ei, 1);
        avg_e_len += (V.row(v1) - V.row(v2)).norm();
    }
    avg_e_len /= static_cast<double>(EV.rows());


    // --- create mask coefficients
    VectorXd coefficients; coefficients.setOnes(nV); // set dummy values to coefficients
    if (with_smoothened_boundary) { // then (fake) boundary vertices are completely fixed, all other vertices move freely
        for (int vi=0; vi<V.rows(); ++vi)
            if (boundary_vis[vi])
                coefficients[vi] = 0;
    }
    else // then create smooth 'mask' coefficients where vertices close to the boundary move less and less.
    {
        if (boundaryPolyline.rows() > 0)
        {
            VectorXd bndry_dist(nV);
            for (int vi = 0; vi < nV; ++vi)
            {
                Vector3d pt = V.row(vi);
                bndry_dist(vi) = (MeshHelpers::projectToPolyline(pt, boundaryPolyline) - pt).norm();
            }
            double threshold = 3.5 * avg_e_len;
            for (int vi = 0; vi < nV; ++vi) {
                if (bndry_dist(vi) < threshold){
                    coefficients[vi] = Helpers::remap_unbound(bndry_dist(vi), 0.0, threshold, 0.0, 1.0); // coefficients go from 1.0 to 0.0 as we approach the boundary
                    coefficients[vi] *= coefficients[vi]; // squared
                }
                else{
                    coefficients[vi] = 1.0;
                }
                    
            }
        }
    }
	// Helpers::write_vector_to_txt(coefficients, DATA_PATH + np3dp->output_folder + string("coefficients.txt"));


    // --- create adjacency lists
    vector<vector<int>> Adjacency(nV);
    vector<vector<int>> Boundary_adjacency(nV);

    auto edge_is_on_boundary = [](int ei, const MatrixXi& EF)->bool { if (EF(ei, 0) == -1 || EF(ei, 1) == -1) return true; return false; };

    for (int ei = EV.rows()-1; ei >=0 ; --ei) {
	        int v1 = EV(ei, 0);
	        int v2 = EV(ei, 1);
        if (!edge_is_on_boundary(ei, EF)) { 
    		Adjacency[v1].push_back(v2);
    		Adjacency[v2].push_back(v1);
        }
        else {
    		Boundary_adjacency[v1].push_back(v2);
    		Boundary_adjacency[v2].push_back(v1);
        }
    }

    // --- debug --->
	//Helpers::write_VV_to_txt(DATA_PATH + np3dp->output_folder + "VV.txt", Adjacency);
	//Helpers::write_VV_to_txt(DATA_PATH + np3dp->output_folder + "VV_boundary.txt", Boundary_adjacency);
    // <---


    // --- construct uniform Laplacians (using adjacency lists) : nVxnV
    SparseMatrix<double> L;
    {
		vector<Triplet<double>> L_triplets;
		for (int vi = 0; vi < nV; ++vi)
		{
		    // off diagonal
		    for (int nv : Adjacency[vi])
		        L_triplets.emplace_back(vi, nv, 1.0 / static_cast<double>(Adjacency[vi].size()));
		    // on the diagonal
		    L_triplets.emplace_back(vi, vi, -1.0);
		}
		L.conservativeResize(nV, nV);
		L.setFromTriplets(L_triplets.begin(), L_triplets.end());
    }

	SparseMatrix<double> L_boundary;
    {
		vector<Triplet<double>> L_triplets;
		for (int vi = 0; vi < nV; ++vi)
		{
            if (!Boundary_adjacency[vi].empty())
            {
			    // off diagonal
			    for (int nv : Boundary_adjacency[vi])
			    {
			        L_triplets.emplace_back(vi, nv, 1.0 / static_cast<double>(Boundary_adjacency[vi].size()));
			    }
				// on the diagonal
				L_triplets.emplace_back(vi, vi, -1.0);
            }
		}
		L_boundary.conservativeResize(nV, nV);
		L_boundary.setFromTriplets(L_triplets.begin(), L_triplets.end());
    }

    
    // --- get per-vertex mass matrix (VxV)
    VectorXd A; // per face areas of triangles
    igl::doublearea(V, T, A);
    A *= 0.5;
    A = MeshHelpers::FtoV(V, T) * A; // convert to per-vertex areas
    SparseMatrix<double> M;
    Helpers::diag(A, M);

    SparseMatrix<double> Id(nV, nV); Id.setIdentity();

    // --- get delta step size find the smallest non-zero eigenvalue of the system L* X = s * M * X 
    double delta = 0.01; // default value in case the eigenvalue problem doesnt converge
    double boundary_delta = 0.5;
	{
        // make positive definite
    	L = L * -1.0; 
        L_boundary = L_boundary * -1.0;

        // --- find smallest positive eigenvalue using matlab
        const auto get_step_size = [](const SparseMatrix<double>& M, const SparseMatrix<double>& L, double default_value)->double {

            Engine* engine;

            MatrixXd EV; // Eigenvectors of the laplacian (w. mass matrix)
            MatrixXd ED; // Eigenvalues of the laplacian (w. mass matrix)

            // Launch MATLAB
            igl::matlab::mlinit(&engine);

            // Send Laplacian matrix to matlab
            igl::matlab::mlsetmatrix(&engine, "L", L);

            // Send mass matrix to matlab
            igl::matlab::mlsetmatrix(&engine, "M", M);

            // Extract the first 5 eigenvectors
            igl::matlab::mleval(&engine, "[EV,ED] = eigs(L,M,5,'smallestabs')");
            // Turn eigenvalue diagonal matrix into a vector
            igl::matlab::mleval(&engine, "ED=diag(ED)");

            // Retrieve the result
            igl::matlab::mlgetmatrix(&engine, "EV", EV);
            igl::matlab::mlgetmatrix(&engine, "ED", ED);

            for (int i = 0; i < ED.rows(); ++i)
            {
                if (ED(i, 0) > 1e-6)
                {
                    return 0.001 / ED(i, 0);
                }
            }
            return default_value; // just return default value
        };

        delta = get_step_size(M, L, delta);
        // boundary_delta = get_step_size(Id, L_boundary, boundary_delta);
        cout << "Step sizes: delta = " << delta << ", boundary_delta = " << boundary_delta << endl;

        // return to original sign
        L = L * -1.0; 
        L_boundary = L_boundary * -1.0;
	}
    

    // --- smoothening iterations
    cout << "--- Smoothing quad mesh with " << iterations << " iterations, step size = " << delta << ", and projecting it back on original surface." << endl;

	//// --- explicit steps
	//for (int iteration = 0; iteration < iterations; ++iteration) {
	//   V = (Id + delta * coefficients.asDiagonal() * L)* V;
	//}

    // --- implicit steps
	SparseLU<SparseMatrix<double>> solver(M - delta * coefficients.asDiagonal() * L); //Id - delta * coefficients.asDiagonal() * L
    if (solver.info() != Eigen::Success) { cerr << "Could not compute implicit solver. Exiting without smoothing." << endl; return; }

    SparseLU<SparseMatrix<double>> solver_boundary;
	if (with_smoothened_boundary){
        solver_boundary.compute(Id - boundary_delta * L_boundary); //Id - delta * coefficients.asDiagonal() * L
        if (solver_boundary.info() != Eigen::Success) { cerr << "Could not compute implicit solver. Exiting without smoothing." << endl; return; }
    }

	for (int iteration = 0; iteration < iterations; ++iteration) {
        cout << "Smoothening iteration " << iteration << endl;

        if (with_smoothened_boundary) // first do laplacian smoothing on the boundary only
        {
			for (int i=0; i<5; ++i)
			{
                // smoothen
                V = solver_boundary.solve(V).eval();

                // project to original boundary
				for (int vi = 0; vi < nV; ++vi) 
					if (boundary_vis[vi] == 1)
						V.row(vi) = MeshHelpers::projectToPolyline(V.row(vi), boundaryPolyline);
			}
            //MeshHelpers::save_mesh_data(V, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);
        }

	    V = solver.solve(M * V).eval();
		//MeshHelpers::save_mesh_data(V, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);

	    // --- project to original mesh surface + boundary
	    VectorXd sqrD; // #P list of smallest squared distances
	    VectorXi I; // #P list of primitive indices corresponding to smallest distances
	    MatrixXd C; // #P by 3 list of closest points
	    //igl::point_mesh_squared_distance(V, Voriginal, T, sqrD, I, C);
	    igl::point_mesh_squared_distance(V, Triangle_V, Triangle_F, sqrD, I, C);
	    V = C;
        //MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);


	    // --- project boundary to original boundary
        {
            for (int vi = 0; vi < nV; ++vi)
                if (boundary_vis[vi] == 1)
                    V.row(vi) = MeshHelpers::projectToPolyline(V.row(vi), boundaryPolyline);
        }


        //MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);
    }
}


///////////////////////////////////
/// --- Rewiring 
///////////////////////////////////


void QuadMesh::rewire_strip(VectorXi& strip, const array<vector<int>, 2>& strip_vis, bool strip_is_closed, bool rewire_forward)
{
    if (!has_updated_mesh_data()) {cerr << "QuadMeshEditing::rewire_strip() requires updated mesh data" << endl; throw;}
    if (strip_is_closed) if (strip_vis[0].size() != strip_vis[1].size()) {cerr << "Attention! Selected strip is closed, but visA.size() != visB.size()" << endl; throw;}

    std::cout << "Rewiring selected strip... ";

    // remember previous data (enables Undo operation)
    remember_previous_data();

    // saving original quad mesh data (for debugging)
    // MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);

    int nF = F.rows();

    // --- decide on orientation if there are selected_vis to be aligned
    if (selection.selected_vis.size() == 2) 
    {
        if (selection.selected_vis[0] != selection.selected_vis[1]){ // TODO: figure out automatically rewiring direction also when selected_vis[0] == selected_vis[1]
	        rewire_forward = decide_rewiring_direction(strip_vis[0], strip_vis[1]);
	        cout << "Rewiring direction : "; if (rewire_forward) cout << "Forward" << endl; else cout << "Backwards" << endl;
        }
    }

    const vector<int>& visA = rewire_forward ? strip_vis[0] : strip_vis[1];
    const vector<int>& visB = rewire_forward ? strip_vis[1] : strip_vis[0];

    // --- find alignment of two groups of vertices
    int i = 0;
    int j = -1;
    int fi_common = -1;
    const bool found = find_alignment_of_two_vertex_groups(strip, visA, visB, i, j, fi_common);
    if (!found) { cerr << "Could not find any connection between visA and visB. Cannot rewire strip" << endl; return; }

    // update i, j so that i==0 and j is updated accordingly to maintain alignment 
    ++j; if (j == visB.size()) j = 0; // return to 0
    while (i>0){
	    --i;
        --j; if (j == 0) break;
    }


    // --- remember order in which vertices should be added on the faces
    bool order1 = false;
    {
        const Vector3d cross1 = Vector3d(V.row(visB[j]) - V.row(visA[i])).cross(Vector3d(V.row(visA[i+1]) - V.row(visA[i])));
        const Vector3d cross2 = Vector3d(V.row(F(fi_common, 2)) - V.row(F(fi_common, 1))).cross(Vector3d(V.row(F(fi_common, 0)) - V.row(F(fi_common, 1))));
        const double dir = cross1.dot(cross2);
        if (dir > 0)
            order1 = true;
    }


    // --- remove strip faces
    for (int fi = nF - 1; fi >= 0; --fi) {
        if (strip[fi]) {
            Helpers::removeRow(F, fi);
            Helpers::removeRow(D, fi);
        }
    }


    // --- add new faces
    // first triangle
    if (!strip_is_closed){
        if (j > 0){
            int fi = F.rows();
            VectorXi face; face.setConstant(F.cols(), -1);

            F.conservativeResize(fi + 1, F.cols());
            D.conservativeResize(fi + 1);

            if (order1) { face[0] = visA[i]; face[1] = visB[j-1]; face[2] = visB[j];}
            else        { face[0] = visA[i]; face[1] = visB[j];   face[2] = visB[j-1];}

            F.row(fi) = face;
            D[fi] = 3;
        }
    }

    
    while (i<visA.size())
    {
	    int fi = F.rows();
        VectorXi face; face.setConstant(F.cols() , -1);

        //cout << "i :" << i << " , j : " << j << endl;
        if (i+1 < visA.size() && j+1 < visB.size() || strip_is_closed) // all intermediate triangles
        {
            F.conservativeResize(fi + 1, F.cols());
            D.conservativeResize(fi + 1);

            if (order1){face[0] = visA[i]; face[1] = visB[j%visB.size()];     face[2] = visB[(j+1)%visB.size()]; face[3] = visA[(i+1)%visA.size()];}
            else       {face[0] = visA[i]; face[1] = visA[(i+1)%visA.size()]; face[2] = visB[(j+1)%visB.size()]; face[3] = visB[j%visB.size()];}
            
	        F.row(fi) = face;
            D[fi] = 4;
            //cout << face << endl;
            ++i; ++j;
        }
        else
        {
            if (!strip_is_closed)
            {
                // last triangle
                if (i+1 < visA.size() && j < visB.size())
                {
                    F.conservativeResize(fi + 1, F.cols());
                    D.conservativeResize(fi + 1);

                    if (order1) {face[0] = visA[i]; face[1] = visB[j];   face[2] = visA[i+1];}
                    else              {face[0] = visA[i]; face[1] = visA[i+1]; face[2] = visB[j];}

                    F.row(fi) = face;
                    D[fi] = 3;
                    //cout << face << endl;
                    ++i; ++j;
                }

                else if (i < visA.size() && j+1 < visB.size())
                {
                    F.conservativeResize(fi+1, F.cols());
                    D.conservativeResize(fi+1);

                    if (order1) {face[0] = visA[i]; face[1] = visB[j];   face[2] = visB[i+1];}
                    else              {face[0] = visA[i]; face[1] = visB[i+1]; face[2] = visB[j];}

                    F.row(fi) = face;
                    D[fi] = 3;
                    //cout << face << endl;
                    ++i; ++j;
                }
            }
            break;
        }
    }
    strip.setZero(F.rows());


    std::cout << "Completed." << endl;

    // ------------------------------------------------------------------------------ Cleaning up --->
    cleanup();

    save_data();

    if (!has_updated_mesh_data()) {cerr << "Attention! QuadMeshEditing::rewire_strip() finished without updated mesh data" << endl; throw;}
}

bool QuadMesh::find_alignment_of_two_vertex_groups(const VectorXi& strip, const vector<int>& visA, const vector<int>& visB, int& i, int&j, int& fi_common)
{

    while (j == -1 && i < visA.size()) {
        int viA = visA[i];

        if (!boundaryVertices[viA]) {
            // find all direct neighbors of viA, and check if any of them are in viB. if so find the index of the neighbor in viB
            for (int fi = 0; fi < F.rows(); ++fi) {
                if (strip[fi]) {
                    for (int k = 0; k < D[fi]; ++k) {
                        if (F(fi, k) == viA) {
                            int nvi_prev = F(fi, (D[fi] + k - 1) % D[fi]);
                            int nvi_next = F(fi, (k + 1) % D[fi]);
                            if (!boundaryVertices[nvi_prev])
                                if (Helpers::in_vector(visB, nvi_prev)) { j = Helpers::element_index(visB, nvi_prev); fi_common = fi; return true; }
                            if (!boundaryVertices[nvi_next])
                                if (Helpers::in_vector(visB, nvi_next)) { j = Helpers::element_index(visB, nvi_next); fi_common = fi; return true; }
                        }
                    }
                }
            }
        }
        ++i;
    }
    return false;
}

void QuadMesh::get_shortest_path_to_vis(int source, vector<int> target_vis, vector<int> other_vis, const vector<vector<int>>& VV, vector<int>& path, bool& path_contains_strip)
{
    set<int> targets; for (int vi : target_vis) targets.insert(vi);
    bool found = QuadMesh::shortest_path(source, targets, VV, path);
    if (!found) throw invalid_argument("Attention! Could not find shortest path from source : " + to_string(source) + " to strip vis");

    assert(Helpers::in_vector(target_vis, path[0])); // last point should always belong to visA
    int second_last_vi       = path.size() > 1 ? path[1] : -1;
    path_contains_strip = path.size() > 1 ? Helpers::in_vector(other_vis, second_last_vi) : false;
}

bool QuadMesh::decide_rewiring_direction(const vector<int>& visA, const vector<int>& visB)
{
    const vector<vector<int>> VV = adjacency_list_with_Emap(0); // get adjacency list of quad mesh

    int source0 = selection.selected_vis[0];
    int source1 = selection.selected_vis[1];
    assert(source0 != source1);

    set<int> targetsA; for (int vi : visA) targetsA.insert(vi);
    set<int> targetsB; for (int vi : visB) targetsB.insert(vi);

    vector<int> path_A0, path_B0, path_A1, path_B1;
    bool A0_contains_strip, B0_contains_strip, A1_contains_strip, B1_contains_strip;

    get_shortest_path_to_vis(source0, visA, visB, VV, path_A0, A0_contains_strip);
    get_shortest_path_to_vis(source0, visB, visA, VV, path_B0, B0_contains_strip);
    get_shortest_path_to_vis(source1, visA, visB, VV, path_A1, A1_contains_strip);
    get_shortest_path_to_vis(source1, visB, visA, VV, path_B1, B1_contains_strip);

    int vA_closest = -1, vB_closest = -1;

    if      (A0_contains_strip && !B0_contains_strip) { vB_closest = path_B0[0]; vA_closest = path_A1[0]; }
    else if (!A0_contains_strip && B0_contains_strip) { vA_closest = path_A0[0]; vB_closest = path_B1[0]; }

    if (vA_closest == -1 || vB_closest == -1)
    {
        if      (A1_contains_strip && !B1_contains_strip){ vB_closest = path_B1[0]; vA_closest = path_A0[0]; }
		else if (!A1_contains_strip && B1_contains_strip){ vA_closest = path_A1[0]; vB_closest = path_B0[0];}
    }

    if (vA_closest == -1 || vB_closest == -1) // then check distances
    {
        if (path_A0.size() > path_B0.size()){
            if (path_A1.size() > path_B1.size()) {cerr << "In strip rewiring: dist_A_0 > dist_B_0 AND dist_A_1 > dist_B_1. This should never happen" << endl; throw;}
            vB_closest = path_B0[0];
            vA_closest = path_A1[0];
        }
        else{
            vA_closest = path_A0[0];
            vB_closest = path_B1[0];
        }
    }

    int i_closest = Helpers::element_index(visA, vA_closest); // Attention! It can be that due to the one side having more boundary vertices, this is not always 100% correct for giving orientation! You need to check for alignment between visA and visB
    int j_closest = Helpers::element_index(visB, vB_closest);

    if (i_closest < j_closest)
        return true;
    else
        return false;
}


///////////////////////////////////
/// --- Subdivide 
///////////////////////////////////

namespace SubdivisionHelpers {

    VectorXi get_strip_edges_with_correct_color(const vector<int>& strip_faces, const VectorXi& D, const VectorXi& Emap, const MatrixXi& FE, int color) {
        VectorXi E; E.setZero(Emap.size());

        for (int fi : strip_faces) {
            vector<int> eis;
            for (int k = 0; k < D[fi]; ++k) {
                int ei = FE(fi, k);
                if (Emap[ei] == color)
                    eis.push_back(ei);
            }
            assert(!eis.empty());
            if (eis.size() == 1) { //  && D[fi] > 3 then also add a boundary edge
                for (int k = 0; k < D[fi]; ++k) {
                    int ei = FE(fi, k);
                    if (Emap[ei] == -1) {
                        if (D[fi] == 3) // on triangle face just add whichever edge is on the boundary
                        {
                            eis.push_back(ei);
                            break;
                        }
                        else // only add edge if the prev and next edges have not been added
                        {
                            if (eis[0] != FE(fi, (k - 1 + D[fi])% D[fi]) && eis[0] !=  FE(fi, (k + 1)% D[fi]))
                            {
                                eis.push_back(ei);
                                break;
                            }
                        }
                    }
                }
            }
            if (eis.size() < 2 && D[fi] > 3) {cerr << "eis.size() != 2 for fi : " + to_string(fi) + " with D[fi] = " + to_string(D[fi]) << endl; throw;}
            // if eis.size()>2, then the 3rd edge is being ignored

            E[eis[0]] = 1;
            if (eis.size() > 1) E[eis[1]] = 1;
        }
        return E;
    }


    void add_new_vertices_on_Esubd(MatrixXi& E_to_Nvi, MatrixXd& V, MatrixXd& UV, const VectorXi& E, const MatrixXi& EV, const vector<int>& eis_start_vi, int N) {
        int nE = EV.rows();
        if (E_to_Nvi.rows() != nE) E_to_Nvi.setConstant(nE, N - 1, -1);

        for (int ei = 0; ei < nE; ++ei) {
            if (E[ei]) {
                if (E_to_Nvi(ei, 0) == -1) { // if no nvis have been added on this ei yet
                    int v1 = EV(ei, 0);
                    int v2 = EV(ei, 1);
                    if (Helpers::in_vector(eis_start_vi, v1))
                        std::swap(v1, v2);

                    int nvi = V.rows();
                    V.conservativeResize(nvi + N - 1, 3);
                    UV.conservativeResize(nvi + N - 1, 2);
                    // np3dp->FixedVis.conservativeResize(nvi + N - 1);

                    double step = 1.0 / static_cast<double>(N);
                    for (int k = 0; k < N - 1; ++k) {
                        double t = (k + 1) * step;
                        V.row(nvi + k) = t * V.row(v1) + (1 - t) * V.row(v2);
                        UV.row(nvi + k) = t * UV.row(v1) + (1 - t) * UV.row(v2);
                        // np3dp->FixedVis[nvi+k] = 0;
                        E_to_Nvi(ei, k) = nvi + k;
                    }
                }
            }
        }
    }


    void divide_D4_face(int fi, MatrixXi& F, VectorXi& D, int N, int nfi, const MatrixXd& V, const VectorXi& old_fi_vis,
        int D_fi, const MatrixXi& E_to_Nvi, const VectorXi& E, const MatrixXi& FE, const MatrixXi& EV) { // split in N quads

		// --- find current face data
        int k = 0;
        int ei = -1, prev_ei = -1, next_ei = -1;
        int v0 = -1, v1 = -1, v2 = -1, v3 = -1;
        for (k = 0; k < D_fi; ++k) {
            ei = FE(fi, k);
            if (!E[ei]) {

                // prev_ei and next_ei are the edges to be subdivided
                prev_ei = FE(fi, (D_fi + k - 1) % D_fi);
                next_ei = FE(fi, (k + 1) % D_fi);
                if (!E[prev_ei]) {
                    // Helpers::write_matrix_to_txt(F, output_folder + "current_F.txt");
                    // Helpers::write_matrix_to_txt(V, output_folder + "current_V.txt");
	                cerr << "E[prev_ei] == 0, this shouldn't happen! prev_ei = " + to_string(prev_ei) + " , vis : " + to_string(EV(prev_ei, 0)) + " , " + to_string(EV(prev_ei, 1))<< ". Ei vis : " << EV(ei, 0) << " , " << EV(ei, 1) <<endl; throw;
                }
                if (!E[next_ei]){
                    // Helpers::write_matrix_to_txt(F, output_folder + "current_F.txt");
                    // Helpers::write_matrix_to_txt(V, output_folder + "current_V.txt");
	                cerr << "E[next_ei] == 0, this shouldn't happen! next_ei = " + to_string(prev_ei) + " , vis : " + to_string(EV(next_ei, 0)) + " , " + to_string(EV(next_ei, 1)) << ". Ei vis : " << EV(ei, 0) << " , " << EV(ei, 1) << endl; throw;
                }

                v0 = old_fi_vis[k];
                v1 = old_fi_vis[(k + 1) % D_fi];
                v2 = old_fi_vis[(k + 2) % D_fi];
                v3 = old_fi_vis[(k + 3) % D_fi];

                if ((V.row(v0) - V.row(E_to_Nvi(prev_ei, 0))).squaredNorm() <= (V.row(v3) - V.row(E_to_Nvi(prev_ei, 0))).squaredNorm())
                    break;
            }
        }
        if (ei == -1) {cerr << "No edge in face fi = " + to_string(fi) + " with D = 4 has an edge with !E[ei]" << endl; throw;}

        // --- fill in new faces
        RowVectorXi nf; nf.setConstant(F.cols(), -1);
        for (int j = 0; j < N; ++j)
        {
            // first (to second last) face                        // second (to last) face
            if (j == 0)  nf(0) = v0;                             else nf(0) = E_to_Nvi(prev_ei, j - 1);
            if (j == 0)  nf(1) = v1;                             else nf(1) = E_to_Nvi(next_ei, j - 1);
            if (j < N - 1) nf(2) = E_to_Nvi(next_ei, j);    else nf(2) = v2;
            if (j < N - 1) nf(3) = E_to_Nvi(prev_ei, j);    else nf(3) = v3;

            if (nf[0] == -1 || nf[1] == -1 || nf[2] == -1 || nf[3] == -1){
                // Helpers::write_matrix_to_txt(F, output_folder + "current_F.txt");
                // Helpers::write_matrix_to_txt(V, output_folder + "current_V.txt");
	            cerr << "While creating new face fi = " << nfi + j << " with D[fi] = 4, there is a vi = -1. New face vertices = " << nf << endl;
            	throw;
            }
            F.row(nfi + j) = nf;
            D[nfi + j] = 4;
        }
    }


    void divide_D5_face(int fi, MatrixXi& F, VectorXi& D, int N, int nfi, const VectorXi& old_fi_vis, int D_fi, const MatrixXi& E_to_Nvi, const VectorXi& E, const MatrixXi& FE) {
        // split in N-1 quads and one D=5 face

		// --- find current face data
        int k;
        int ei = -1; // find ei that has subdivision and has no vertices on the boundary
        int prev_ei = -1, next_ei = -1; // prev_ei and next_ei are the edges to be subdivided
        bool found = false;
        for (k = 0; k < D_fi; ++k) {
            ei = FE(fi, k);
            prev_ei = FE(fi, (D_fi + k - 1) % D_fi);
            next_ei = FE(fi, (k + 1) % D_fi);
            if (!E[ei] && E[prev_ei] && E[next_ei]){
                found = true;
                break;
            }
        }
        if (!found) { cerr << "No edge in face fi = " + to_string(fi) + " with D = 5 has an edge with a subdivision nvi and no vertex on the boundary" << endl; throw; }

    	// --- fill in new faces
        RowVectorXi nf; nf.setConstant(F.cols(), -1); // stores new faces
        for (int j = 0; j < N - 1; ++j) // this part is the same as for D=4
        {
            // first (to second last) face                               // second (to last) face
            if (j == 0)  nf(0) = old_fi_vis[k];                     else nf(0) = E_to_Nvi(prev_ei, j - 1);
            if (j == 0)  nf(1) = old_fi_vis[(k + 1) % D_fi];        else nf(1) = E_to_Nvi(next_ei, j - 1);
            if (j < N - 1) nf(2) = E_to_Nvi(next_ei, j);    else nf(2) = old_fi_vis[(k + 2) % D_fi];
            if (j < N - 1) nf(3) = E_to_Nvi(prev_ei, j);    else nf(3) = old_fi_vis[(k + 3) % D_fi];

            F.row(nfi + j) = nf;
            D[nfi + j] = 4;
        }

        // add last face as D=5 (the only real difference from the D4 function)
        nf.setConstant(F.cols(), -1); // reset
        nf(0) = E_to_Nvi(prev_ei, N - 2);
        nf(1) = E_to_Nvi(next_ei, N - 2);
        nf(2) = old_fi_vis[(k + 2) % D_fi];
        nf(3) = old_fi_vis[(k + 3) % D_fi];
        nf(4) = old_fi_vis[(k + 4) % D_fi];
        for (int i = 5; i < D_fi; ++i) // only difference for faces that have D_fi > 5 is here
        {
            nf(i) = old_fi_vis[(k + i) % D_fi];
        }
        F.row(nfi + N - 1) = nf;
        D[nfi + N - 1] = D_fi;
    }


    void divide_D3_face(int fi, MatrixXi& F, VectorXi& D, int N, int nfi, const VectorXi& old_fi_vis,
        int D_fi, const MatrixXi& E_to_Nvi, const VectorXi& E, const MatrixXi& FE) {
        // split in N-1 quads and one triangle face
        int k;
        int ei = -1, next_ei = -1;
        for (k = 0; k < D_fi; ++k) {
            ei = FE(fi, k);
            next_ei = FE(fi, (k + 1) % D_fi);
            if (E[ei] && E[next_ei]) // if we are on an edge with E[ei] AND its next edge is also subdivided (either as boundary, or as part of the strip)
                break;
        }
        if (ei == -1) throw invalid_argument("No edge in face fi = " + to_string(fi) + " with D = 3 has an edge with E[ei] == 1");

        const int ei_v1 = old_fi_vis[k]; // v1 on ei
        const int ei_v2 = old_fi_vis[(k + 1) % D_fi]; // v2 on ei
        const int vi_other = old_fi_vis[(k + 2) % D_fi];

        RowVectorXi nf; nf.setConstant(F.cols(), -1); // stores new faces
        for (int j = 0; j < N-1; ++j) // this part is the same as for D=4
        {
            // first (to second last) face                        // second (to last) face
            if (j == 0)  nf(0) = ei_v1;                     else nf(0) = E_to_Nvi(ei, j-1);
                         nf(1) = E_to_Nvi(ei, j);
                         nf(2) = E_to_Nvi(next_ei, j);
        	if (j == 0)  nf(3) = vi_other;                  else nf(3) = E_to_Nvi(next_ei, j-1);

            //cout << nf << endl;
        	F.row(nfi + j) = nf;
            D[nfi + j] = 4;
        }

        // add last triangle
        nf.setConstant(F.cols(), -1); // reset
        nf(0) = E_to_Nvi(ei, N - 2);
        nf(1) = ei_v2;
        nf(2) = E_to_Nvi(next_ei, N - 2);
        //cout << nf << endl;
        F.row(nfi + N - 1) = nf;
        D[nfi + N - 1] = 3;
    }


    bool is_elongated_strip(const RowVectorXi& strip, const VectorXd& A, double avg_area)
    {
        // get average area of each face on current strip
        double current_avg_area = 0;
        for (int fi=0; fi<strip.size(); ++fi) {
            if (strip[fi]){
                // if (A[fi] > 8 * avg_area) {// if the area of any face is more than 3 * avg_area -> subdivide
                //     cout << "---  SUBD REASON : single face too large, fi = " << fi << endl;
                // 	return true;
                // }
	            current_avg_area += A[fi];
            }
        }
        current_avg_area /= static_cast<double>(strip.sum());

        // if the average area per quad is more than 2 * avg_area -> subdivide
        if (current_avg_area > 2.0 * avg_area)
        {
	        cout << "--- SUBD REASON : overall area too large. current_avg_area/avg_area = " << current_avg_area/avg_area << endl;
            return true;
        }

        return false;
    }

}


void QuadMesh::subdivide_strip(const VectorXi& strip, int direction, bool save_results)
{
    if (!has_updated_mesh_data()) { cerr << "QuadMeshEditing::subdivide2_strip() requires updated mesh data" << endl; throw; }

    std::cout << "Subdividing selected strip... ";

    // -------- get data
    remember_previous_data(); // remember previous data (enables Undo operation)

    assert(Emap.rows() == EV.rows());

    int e_remember_v1, e_remember_v2, e_remember_Emap;
    vector<VectorXi*> F_data;

    // -------- subdivide
    subdivide_strip(V, F, D, Emap, UVcoords, F_data, strip, direction, e_remember_v1, e_remember_v2, e_remember_Emap);

    // -------- cleanup
	if(V.rows() != UVcoords.rows()) {cerr << "After subdivision : Vquad.rows() != UVcoordsQuad.rows()" << endl; throw;}

    std::cout << "Completed." << endl;

    // --- update data
    // get back Fdata
    // Fmap = F_data[0];
    // F_connections_map = F_data[1];
    // F_connections_tapering_map = F_data[2];

    update_quad_mesh_data();
    create_quad_Emap();

    // --- make sure new Emap matches old by comparing 
    for (int ei=0; ei<EV.rows(); ++ei){
        int v1 = EV(ei, 0);
        int v2 = EV(ei, 1);
        if (min(v1, v2) == min(e_remember_v1, e_remember_v2) && max(v1, v2) == max(e_remember_v1, e_remember_v2)) // if we are on the same edge
            if (e_remember_Emap != Emap[ei])
                MeshHelpers::flip_Emap(Emap);
    }

    // --- save 
    if (save_results)
        save_data();

    if (!has_updated_mesh_data()) {cerr << "QuadMeshEditing::subdivide2_strip() finished without updating mesh data" << endl; throw;}
}


void QuadMesh::subdivide_strip(MatrixXd& V, MatrixXi& F, VectorXi& D, VectorXi& Emap, MatrixXd& UV, const vector<VectorXi*>& F_based_data,
									  const VectorXi& strip, int direction, int& e_remember_v1, int& e_remember_v2, int& e_remember_Emap)
{
    int nF = F.rows();

    int N = 2; // number of subdivisions (i.e. each face is split in N faces)

    int other_direction = static_cast<int>(!static_cast<bool>(direction));

    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);


    // --- find edges with the correct color for every face that belongs in the strip
    vector<int> strip_faces; for (int fi = 0; fi < nF; ++fi) if (strip[fi]) strip_faces.push_back(fi);
    VectorXi E = SubdivisionHelpers::get_strip_edges_with_correct_color(strip_faces, D, Emap, FE, other_direction);


    // --- add one new vi for every ei with E[ei] == 1
    MatrixXi E_to_Nvi; //stores which new vertices where added in which edge 
    SubdivisionHelpers::add_new_vertices_on_Esubd(E_to_Nvi, V, UV, E, EV, vector<int>(), N);
	//Helpers::write_matrix_to_txt(V, DATA_PATH + np3dp->output_folder + "V_subd.txt");

    // --- find ei that is not in strip and save its Emap + vis
    int ei=0;
    while(true){
	    int f1 = EF(ei, 0);
        int f2 = EF(ei, 1);
        if (f1 >=0 && f2 >= 0)
			if (!strip[f1] && !strip[f2] && Emap[ei] != -1)
				break;
    	++ei;
    }

    e_remember_v1 = EV(ei, 0);
    e_remember_v2 = EV(ei, 1);
    e_remember_Emap = Emap[ei];


    // --- Replace each strip face with two new faces
    for (int fi = nF - 1; fi >= 0; --fi) {
        if (strip[fi]) {

            // --- remember fi data before they are removed from F and D
            VectorXi old_fi_vis = F.row(fi);
            const int D_fi = D[fi];

            // --- remove fi from Face data
            Helpers::removeRow(F, fi);
            Helpers::removeRow(D, fi);

            vector<int> remember_vals;
            for (VectorXi* F_data : F_based_data){
                remember_vals.push_back((*F_data)[fi]);
                Helpers::removeRow(*F_data, fi);
            }

            // --- create N new faces
            int nfi = F.rows();
            F.conservativeResize(nfi + N, F.cols());
            D.conservativeResize(nfi + N);

            for (int jj = 0; jj < F_based_data.size(); ++jj) {
                VectorXi* F_data = F_based_data[jj];
                if (F_data->rows() == F.rows() - N) {
                    F_data->conservativeResize(nfi + N);
                    for (int i = 0; i < N; ++i) {
                        (*F_data)[nfi + i] = remember_vals[jj];
                    }
                }
            }

            if (D_fi == 4) { // then split in N quads
                SubdivisionHelpers::divide_D4_face(fi, F, D, N, nfi, V, old_fi_vis, D_fi, E_to_Nvi, E, FE, EV);
            }
            else if (D_fi >= 5) { // then split in one D=5 and one quad 
                SubdivisionHelpers::divide_D5_face(fi, F, D, N, nfi, old_fi_vis, D_fi, E_to_Nvi, E, FE);
            }
            else if (D_fi == 3) {// then split in two triangles
                SubdivisionHelpers::divide_D3_face(fi, F, D, N, nfi, old_fi_vis, D_fi, E_to_Nvi, E, FE);
            }
            else{
	            cerr << "I don't know how to subdivide face : " << fi << ", with D[fi] = " << D_fi << " and vis : " << old_fi_vis << endl; throw;
            }
        }
    }
}



void QuadMesh::subdivide_whole_quad_mesh(bool only_elongated_strips)
{
    // remember previous data (enables Undo operation)
    remember_previous_data();

    int side = 0; // it doesn't matter which side we consider

    cleanup();
    // save_data();

    // find average quad size
    VectorXd A; double avg_area = 0;
    if (only_elongated_strips){
	    A = face_areas(V, F, T, TF);
        avg_area = A.mean();
    }

    for (int dir = 0; dir < 2; ++dir)
    {
		// get all strips
		MatrixXi strips; // Strips x F (with values strips(si, fi) = 1/0 if fi belongs/doesn't belong to strip si)
	    VectorXi strips_directions; // Strips x 1 (with values 0/1 for strip in the dominant/subdominant direction)
	    vector<array<vector<int>, 2>> strips_vis;
	    vector<bool> strips_are_closed;
	    vector<vector<int>> StoEis, StoFis; // Strip to Edge indices. For each strip: a vector of ordered edge indices
		MatrixXi FtoS; // Fx2: face to strip (for each strip-direction)
	    vector<vector<set<int>>> VtoS; // Vx2: vertex to strip (for each strip-direction):  VtoS[vi][strip_direction][si]

	    Strips::get_all_strips(V, F, D, boundaryFaces, Emap, strips, strips_directions,strips_vis, strips_are_closed, StoEis, StoFis, FtoS, VtoS, {dir});
		 //Helpers::write_matrix_to_txt(strips, DATA_PATH + np3dp->output_folder + "strips.txt");
	     //Helpers::write_vector_to_txt(strips_directions, DATA_PATH + np3dp->output_folder + "strips_directions.txt");

		int nS = strips.rows();

	    for (int si=0; si<nS; ++si)
	    {
            int nF = F.rows();

            bool subdivide = only_elongated_strips ? SubdivisionHelpers::is_elongated_strip(strips.row(si), A, avg_area) : true;

            if (subdivide)
            {
		        // subdivide
		        subdivide_strip(strips.row(si), dir, false);
	            int new_nF = F.rows();

		        // update strips matrix by removing cols that correspond to the deleted fis ...
		        for (int fi=nF-1; fi >=0; --fi){
	                if (strips(si, fi) == 1) { // if this face was removed
	                    Helpers::removeCol(strips, fi);
	                }
		        }

	            // ... and adding cols that correspond to the new fis.
	            int k = strips.cols();
	            strips.conservativeResize(nS, new_nF);

	            // Set all new cols to 0
	            strips.block(0, k, nS, new_nF - k).setZero();

	            // set new cols of subdivided strip to 1
	            strips.row(si).tail(new_nF - k).setOnes();

            	//Helpers::write_matrix_to_txt(strips, DATA_PATH + np3dp->output_folder + "strips.txt");

                // update face areas
                if (only_elongated_strips){
				    A = face_areas(V, F, T, TF);
			    }

                //MeshHelpers::save_mesh_data(np3dp->Vquad, np3dp->Fquad, np3dp->D, np3dp->Emap, np3dp->EVquad, DATA_PATH + np3dp->output_folder);
                //cout << endl;
            }
	    }
        cleanup();
    }

    save_data();
}



////////////////////////
/// --- (de)serialize and invalidate
////////////////////////

void QuadMesh::serialize(const std::string& serialize_folder) const
{
	igl::serialize(V, DATA_PATH + serialize_folder + "quad_mesh.V.igl");
    igl::serialize(F, DATA_PATH + serialize_folder + "quad_mesh.F.igl");
    igl::serialize(D, DATA_PATH + serialize_folder + "quad_mesh.D.igl");
    igl::serialize(Emap, DATA_PATH + serialize_folder + "quad_mesh.Emap.igl");
    igl::serialize(UVcoords, DATA_PATH + serialize_folder + "quad_mesh.UVcoords.igl");
    igl::serialize(ERibMap, DATA_PATH + serialize_folder + "quad_mesh.ERibMap.igl");

    igl::serialize(ribs_dist, DATA_PATH + serialize_folder + "quad_mesh.ribs_dist.igl");
    igl::serialize(ribs_f0, DATA_PATH + serialize_folder + "quad_mesh.ribs_f0.igl");
    igl::serialize(vis_to_check_for_alignments, DATA_PATH + serialize_folder + "quad_mesh.vis_to_check_for_alignments.igl");
}

void QuadMesh::deserialize(const std::string& serialize_folder)
{
	igl::deserialize(V, DATA_PATH + serialize_folder + "quad_mesh.V.igl");
    igl::deserialize(F, DATA_PATH + serialize_folder + "quad_mesh.F.igl");
    igl::deserialize(D, DATA_PATH + serialize_folder + "quad_mesh.D.igl");
    igl::deserialize(Emap, DATA_PATH + serialize_folder + "quad_mesh.Emap.igl");
    igl::deserialize(UVcoords, DATA_PATH + serialize_folder + "quad_mesh.UVcoords.igl");
    igl::deserialize(ERibMap, DATA_PATH + serialize_folder + "quad_mesh.ERibMap.igl");

    igl::deserialize(ribs_dist, DATA_PATH + serialize_folder + "quad_mesh.ribs_dist.igl");
    igl::deserialize(ribs_f0, DATA_PATH + serialize_folder + "quad_mesh.ribs_f0.igl");
    igl::deserialize(vis_to_check_for_alignments, DATA_PATH + serialize_folder + "quad_mesh.vis_to_check_for_alignments.igl");

    update_quad_mesh_data();
    save_data();
}
 