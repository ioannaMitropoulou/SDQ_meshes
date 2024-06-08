#include "paths_tracing.h"

#include <boost/graph/connected_components.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <directional/polygonal_edge_topology.h>
#include <igl/avg_edge_length.h>
#include "quad_mesh.h"
#include <hedra/polygonal_write_OFF.h>
#include <hedra/triangulate_mesh.h>
#include "strips.h"


namespace PathTracingHelpers// some local helpers
{
    bool find_visited_direct_neighbor(int fi, int acceptable_ei_type, const VectorXi& visited, const VectorXi& Emap, const VectorXi& D,
													   const MatrixXi& FE, const MatrixXi& EF, int& nfi, int& ei)
    {
        // acceptable_ei_type = 1 / 0 for faces that are neighbors through a subdominant/dominant edge

        int m = 0; nfi = -1; ei = -1;
        while (m < D[fi]) {
            ei = FE(fi, m);
            if (Emap[ei] == acceptable_ei_type) {
                nfi = EF(ei, 0) == fi ? EF(ei, 1) : EF(ei, 0);
                if (nfi >= 0)
                    if (visited[nfi])
                    	return true;
            }
            ++m;
        }
        nfi = -1; ei = -1; // reset
        return false;
    }


    double is_edge_intersected(double z1, double z2, double isovalue) { // returns intersection parameter (double [0-1] or Nan)
        double s = (isovalue - z1) / (z2 - z1);
        if (s < 0 || s >= 1) // If not intersected, then put NaN
            s = std::numeric_limits<double>::quiet_NaN();
        return s;
    };


    bool is_regular_face(int fj, const VectorXi& D, const VectorXi& Emap, const MatrixXi& FE) // returns true if the face is a) a quad and b) has two dominant edges
    {
        int count = 0;
        int count_subd = 0;
        if (D[fj] != 4)
            return false;

        for (int m = 0; m < D[fj]; ++m)
        {
            if (Emap[FE(fj, m)] == 0)
                ++count;
            if (Emap[FE(fj, m)] == 1)
                ++count_subd;
        }

        bool has_two_dominant_edges = count == 2;
        bool has_two_subdominant_edges = count_subd == 2;
        return has_two_dominant_edges || has_two_subdominant_edges;
    }


    
    bool face_has_rib(int fj, const VectorXi& D, const VectorXi& ERibMap, const MatrixXi& FE) // returns true if the face is a) a quad and b) has two dominant edges
    {
        for (int m = 0; m < D[fj]; ++m)
            if (ERibMap[FE(fj, m)] == 1)
                return true;

        return false;
    }

    bool neighbor_face_has_rib(int fj, const VectorXi& D, const VectorXi& ERibMap, const MatrixXi& FE, const MatrixXi& EF)
    {
	    for (int m = 0; m < D[fj]; ++m){
            int ei = FE(fj, m);
            int nfi = EF(ei, 0) == fj ? EF(ei, 1) : EF(ei, 0);
            if (nfi >= 0){
                if (face_has_rib(nfi, D, ERibMap, FE))
                    return true;
            }
	    }
        return false;
    }


    bool find_regular_face_on_current_row(int& fj, const MatrixXi& F, const VectorXi& D, const VectorXi& Emap, const MatrixXi& EF, const MatrixXi& FE, VectorXi& temp_visited)
    {
        std::deque<int> temp_R = { fj };
        temp_visited.setZero(F.rows());
        temp_visited[fj] = 1;

        do
        {
            fj = temp_R.front();
            temp_R.pop_front();

            for (int m = 0; m < D[fj]; ++m) {
                int ei = FE(fj, m);
                if (Emap[ei] == 1) {
                    int nfj = EF(ei, 0) == fj ? EF(ei, 1) : EF(ei, 0);
                    if (nfj >= 0) {
                        if (!temp_visited[nfj]) {
                            temp_R.push_back(nfj);
                            temp_visited[nfj] = 1;
                        }
                    }
                }
            }
        } while (!is_regular_face(fj, D, Emap, FE) && !temp_R.empty());

    	return is_regular_face(fj, D, Emap, FE);
    }

    bool is_uncut_vertex(int fi, int m, const MatrixXi& FH, const VectorXi& Emap, const VectorXi& prevH, const VectorXi& HE)
    { // i.e. if the vertex belongs to the uncut quad mesh (with integer UV coordinates)
        int he = FH(fi, m);
        if (Emap[HE[he]] >= 0 && Emap[HE[prevH[he]]] >= 0) // if both the edge before and the edge after the vertex have been colored, then it is uncut
            return true;
        return false;
    };


    bool is_unit_edge(int fi, int m, const MatrixXi& EV, const MatrixXi& FH, const VectorXi& Emap, const VectorXi& prevH, const VectorXi& HE, const VectorXi& D)
    { // if both vertices of the edge are uncut, then the edge is unit norm
        if (is_uncut_vertex(fi, m, FH, Emap, prevH, HE) && is_uncut_vertex(fi, (m + 1) % D[fi], FH, Emap, prevH, HE))
            return true;
        return false;
    };
}


void PathsTracing::trace_continuous_paths_on_quad_mesh(Np3dpContext* np3dp, const vector<int>& Piece_IDs) {
    // --- assign number of tracing paths per dominant strip
    cout << endl << "Deciding on number of paths per strip " << endl;
    for (int ID : Piece_IDs) {
        cout << " ID = " << ID << endl;
        Piece& p = np3dp->pieces[ID];
        assign_number_of_paths_per_strip(p, np3dp->default_paths_density, np3dp->mesh_scale_factor, np3dp->desired_avg_layer_height);

        Helpers::write_vector_to_txt(p.nP_per_strip, Helpers::get_filename(DATA_PATH + np3dp->output_folder, "nP_per_strip", ID, "txt"));
    }


    auto start = chrono::steady_clock::now();

    // --- trace paths
    for (int ID : Piece_IDs) {
        cout << endl << "Paths tracing on piece with ID = " << ID << endl;
        Piece& p = np3dp->pieces[ID];

        p.pathsCollection = Paths::PathCollection(); // remove any existing values
        trace(p, np3dp->default_paths_density);
    }

    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;


    // --- save paths to file
    for (int ID : Piece_IDs) 
        np3dp->pieces[ID].pathsCollection.write_to_file(Helpers::get_filename(DATA_PATH + np3dp->output_folder, "paths_collection", ID, "txt"));
    
}


void PathsTracing::trace(Piece& p, int density) {

    const MatrixXd& V = p.V; // current piece vertices
    const MatrixXi& F = p.F; // current piece faces
    const VectorXi& D = p.D; // current piece face degrees
    const MatrixXd& N = p.faceNormals; // current piece face normals
    const VectorXi& Emap = p.Emap;
    const VectorXi& ERibMap = p.ERibMap;
    const MatrixXd& UV = p.UV;

    Paths::PathCollection& paths_collection = p.pathsCollection;
    MatrixXd& guiding_dirs = p.GuidingDirs;

    // --- Some notes
	// lateral direct neighbors are faces connected through subdominant edges (i.e. Emap[ei] = 1)
	// forward/backward direct neighbors are faces connected through dominant edges (i.e. Emap[ei] = 0)
	// i index (ex. ei, nfi etc) refers to forward/backward neighboring. j index (ex.ej, nfj etc refers to lateral neighboring)

    MatrixXi T; VectorXi TF;
    hedra::triangulate_mesh(D, F, T, TF);

    // --- get edge topology data
    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
    hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);

    auto are_neighboring_faces = [](int f1, int f2, const MatrixXi& FE)->bool
    {
        vector<int> edges_f1, edges_f2;
	    for (int i=0; i<FE.cols(); ++i)
	    {
            if (FE(f1, i) >= 0)
                edges_f1.push_back(FE(f1, i));
            if (FE(f2, i) >= 0)
                edges_f2.push_back(FE(f2, i));
        }
        return !Helpers::common_elements_in_vectors(edges_f1, edges_f2).empty();
    };

    // --- get halfedge datastructure 
    MatrixXi EH, FH; VectorXi VH, HE, HF, HV, nextH, prevH, twinH;
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

    // --- get d_threshold value
    double avg_edge_len = igl::avg_edge_length(V, T);
    double d_threshold = 0.85 * avg_edge_len / double(density);
    //cout << "d_threshold = " << d_threshold << endl;

    // --- pick a random starting quad
    int f_start = 0;
    while (!PathTracingHelpers::is_regular_face(f_start, D, Emap, FE) || PathTracingHelpers::face_has_rib(f_start, D, ERibMap, FE) || PathTracingHelpers::neighbor_face_has_rib(f_start, D, ERibMap, FE, EF) ) {
        // while face is NOT regular, or it HAS a rib, move to next face
        ++f_start;
        if (f_start == D.size() - 1) break;
    }

    if (guiding_dirs.rows() != F.rows()){
	    p.get_guiding_directions();
    }


    // --- start expansion
    std::deque<int> Q = { f_start };

    VectorXi visited;
    visited.setZero(F.rows());


    while(!Q.empty()) {
        int fi = Q.front();
        Q.pop_front();
        if (visited[fi]) continue;


        // --------------------------- Visit entire row of fi and create paths
        std::deque<int> R = { fi }; // Row queue
        Paths::PathCollection current_row_paths; // create new paths for the current row R

        int current_nP = p.nP_per_strip[fi];// density; // number of paths in current row
        current_row_paths.paths.resize(current_nP);
        std::set<int> faces_in_current_row;
        

        while (!R.empty()) {
            int fj = R.front();
            R.pop_front();
            if (visited[fj]) continue;

            faces_in_current_row.insert(fj);


            if (current_nP > 0)
            {
                bool regular_face = PathTracingHelpers::is_regular_face(fj, D, Emap, FE);
                bool first_in_row = current_row_paths.paths[0].segments.empty();

                // --- try to have the first face of the row be a regular face (so that the paths are set with the correct orientation)
                if (!regular_face && first_in_row) {
                    visited[fj] = 0; // set fj as unvisited
                    VectorXi temp_visited;
                    regular_face = PathTracingHelpers::find_regular_face_on_current_row(fj, F, D, Emap, EF, FE, temp_visited);

                    if (!regular_face) { // if there is no regular face on the row
                        cout << "Could not find any regular face in current row where fj = " << fj << " belongs." << endl;
                        int max_ff = -1, max_segs_num = -1;
                        for (int ff = 0; ff < temp_visited.size(); ++ff) { // find which face in the row has the most segments
                            if (temp_visited[ff] == 1) {
                                auto segs = trace_segments_on_irregular_face(ff, current_nP, V, F, D, Emap, UV, FH, HE, prevH, EV, guiding_dirs.row(ff), avg_edge_len);
                                if (static_cast<int>(segs.size()) > max_segs_num) {
                                    max_segs_num = segs.size();
                                    max_ff = ff;
                                    if (max_segs_num == current_nP) break;
                                }
                            }
                        }
                        assert(max_ff >= 0);
                        fj = max_ff;
                        current_row_paths.paths.resize(max_segs_num); // change the number of paths
                        current_nP = max_segs_num;
                    }

                    fi = fj; // also update fi 
                }
                assert(current_row_paths.paths.size() == current_nP);


                visited[fj] = 1; // set as visited
                faces_in_current_row.insert(fj);


                // --- get segments on current face fi (the resulting faces are already looking in the same direction)
                vector<Paths::Segment> segments;
                if (regular_face) // quad faces not cut by boundaries
                {
                    segments = trace_segments_on_quad_face(fj, current_nP, V, FH, HE, nextH, HV, Emap);
                }
                else if (D[fj] == 4) // exceptionally treat as regular
                {
                    segments = trace_segments_on_quad_face(fj, current_nP, V, FH, HE, nextH, HV, Emap);
                }
                else // non-quad or boundary-cut faces
                {
                    segments = trace_segments_on_irregular_face(fj, current_nP, V, F, D, Emap, UV, FH, HE, prevH, EV, guiding_dirs.row(fj), avg_edge_len);
                }


                // --- order segments and add them to the correct paths
                if (!segments.empty())
                {
                    // --- set segments forward_dir 
                    for (auto& seg : segments)
                        seg.forward_dir = guiding_dirs.row(fj);

                    // --- order segments based on forward_dir
                    if (segments.size() > 1) {
                        if ((segments[1].v0 - segments[0].v0).dot(Vector3d(guiding_dirs.row(fj))) < 0) {
                            std::reverse(segments.begin(), segments.end());
                            for (auto& seg : segments)
                                seg.switch_v0_v1();
                        }
                    }

                    // --- orient segments so that forward_dir.cross(avg_segment_dir) || face Normal
                    //if (!regular_face)
                    {
                        Vector3d avg_dir(0,0,0);
                        for (const auto seg : segments)
                            avg_dir += (seg.v1 - seg.v0).normalized();
                        avg_dir.normalize();
                        if (  avg_dir.cross(Vector3d(guiding_dirs.row(fj))).dot(Vector3d(N.row(fj))) < 0 )
                            for (auto& seg : segments)
                                seg.switch_v0_v1();
                    }
                    


                    // --- calculate smooth normals for each segment
                    {
	                    // find prev and next face orthogonal to the strip
	                    auto get_face_centroid = [](int fi, const VectorXi& D, const MatrixXd& V, const MatrixXi& F)->RowVector3d {
	                        RowVector3d cen(0, 0, 0);
                    		for (int m = 0; m < D[fi]; ++m) 
	                            cen += V.row(F(fi, m));
                    		cen /= static_cast<double>(D[fi]);
	                        return cen;
	                    };

	                    int next_fi = -1; int prev_fi = -1;
	                    Vector3d cen_fj = get_face_centroid(fj, D, V, F);
	                    for (int m = 0; m < D[fj]; ++m) {
	                        int ei = FE(fj, m);
	                        if (Emap[ei] == 0) {
	                            int nfi = EF(ei, 0) == fj ? EF(ei, 1) : EF(ei, 0);
	                            if (nfi >= 0) {
	                                Vector3d cen_nfi = get_face_centroid(nfi, D, V, F);
	                                if ((cen_fj - cen_nfi).dot(guiding_dirs.row(fj)) < 0)
	                                    next_fi = nfi;
	                                if ((cen_fj - cen_nfi).dot(guiding_dirs.row(fj)) > 0)
	                                    prev_fi = nfi;
	                            }
	                        }
	                        if (next_fi >=0 && prev_fi>=0) break;
	                    }


	                    // cosine weights
	                    // x -> from 0 to 1,    cos() -> from 1 to -1,    0.25*cos() -> from 0.25 to -0.25,    0.75 + 0.25*cos() -> from 1.0 to 0.5
	                    auto get_param2 = [](double x)->double { return 0.25 + 0.25 * cos(x * 2 * igl::PI); }; // returns values from 0.0 to 0.5, peaks for x = 0.0 and x = 1.0
	                    auto get_param1 = [](double x)->double { return 0.75 - 0.25 * cos(x * 2 * igl::PI); }; // returns values from 0.5 to 1.0, peaks for x= 0.5

	                    // linear weights
	                    //auto get_param2 = [](double x)->double { return 0.5 * abs(1. - 2.0 * x); }; // returns values from 0.0 to 0.5, peaks for x = 0.0 and for x = 1.0
	                    //auto get_param1 = [](double x)->double { return 0.5 + 0.5 * (1 - abs(1. - 2.0 * x)); };// returns values from 0.5 to 1.0, peaks for x = 0.5

	                    int n = segments.size();
	                    for (int si=0; si<n; ++si)
	                    {
		                    double x = static_cast<double>(2*si + 1) / static_cast<double>(2*n);
	                        //double x = static_cast<double>(si) / static_cast<double>(n-1);
	                        double p1 = get_param1(x);
	                        double p2 = get_param2(x);
	                        //cout << "x : " << x << ", p1 : " << p1 << ", p2 " << p2 << endl;
	                        const RowVector3d& normal_of_current_face = N.row(fj);

	                        if (x < 0.5 && prev_fi >= 0) {
	                            segments[si].normal = p1 * normal_of_current_face + p2 * N.row(prev_fi);
	                        }
	                        else if (x>=0.5 && next_fi >= 0) {
	                            segments[si].normal = p1 * normal_of_current_face + p2 * N.row(next_fi);
	                        }
	                        else {
	                            segments[si].normal = normal_of_current_face;
	                        }
	                        segments[si].normal.normalize();
	                    }
                    }



                    // --- deal with irregular faces with more or less segments than current_nP
                    int position_i = 0; // which path index the first segment corresponds to (not 0 only when the segments of current face are less than the number of paths)
                    if (!regular_face) {
                        if (!first_in_row) { // if there are existing segments to compare with


	                        if (segments.size() > current_nP) // --- remove extra segments
	                        {
	                            assert(QuadMesh::is_boundary_face(fj, twinH, D, FH)); // this can only happen on the boundary faces
	                            int n = segments.size() - current_nP; // number of segments to be removed
	                            for (int i = 0; i < n; ++i) {
	                                // check first and last segment, and decide which one should be removed
                                    bool found = false;
                                    double d0 = 1e10; double d_last = 1e10; // very big numbers 

                                    for (int attempt = 0; attempt < 2; ++attempt) {
                                        if (segments[0].parentEdges[1] == current_row_paths.paths[0].ei_first()) {
                                            d0 = (segments[0].v1 - current_row_paths.paths[0].v_first()).norm();
                                            found = true;
                                        } else if (segments[0].parentEdges[0] == current_row_paths.paths[0].ei_last()) {
                                            d0 = (segments[0].v0 - current_row_paths.paths[0].v_last()).norm();
                                            found = true;
                                        }

                                        if (segments.back().parentEdges[1] == current_row_paths.paths[0].ei_first()) {
                                            d_last = (segments.back().v1 - current_row_paths.paths[0].v_first()).norm();
                                            found = true;
                                        } else if (segments.back().parentEdges[0] == current_row_paths.paths[0].ei_last()) {
                                            d_last = (segments.back().v0 - current_row_paths.paths[0].v_last()).norm();
                                            found = true;
                                        }

                                        //if (!found)
                                        //    for (auto& seg : segments)
                                        //        seg.switch_v0_v1();
                                        //else
                                        //    break;
                                        if (found)
                                            break;
                                    }

                                    if (d0 > d_last)
                                        segments.erase(segments.begin()); // remove first element
                                    else
                                        segments.pop_back(); // remove last element
	                            }
                                assert(segments.size() == current_nP);
                                cout << "Removed " << n << " segments from face fj = " << fj << " to match the row" << endl;
	                        }


                            else if (segments.size() < current_nP)  // for faces that have less that #current_nP segments traced
                            {
                                int options = current_nP - segments.size();

                                // ====== find side (beginning or end)
                                bool found = false;
                                bool at_beginning = false;
                                int side_finding_i = -1;

                                while(!found && !segments.empty()) {
	                                options = current_nP - segments.size();
	                                for (int times = 0; times < 2; ++times) {
	                                    for (side_finding_i = 0; side_finding_i <= options; ++side_finding_i) {
	                                        int pi = side_finding_i;

	                                        // at the beginning of current paths
	                                        if (segments[0].parentEdges[1] == current_row_paths.paths[pi].ei_first()) {
	                                            at_beginning = true; found = true;
	                                            break;
	                                        }

	                                        // at the end
	                                        if (segments[0].parentEdges[0] == current_row_paths.paths[pi].ei_last()) {
	                                            at_beginning = false; found = true;
	                                            break;
	                                        }
	                                    }

	                                    if (found)
	                                        break;
	                                    //else if (times == 0) // only the 1st time
	                                    //    for (auto& seg : segments) // flip segments orientation and try again
	                                    //        seg.switch_v0_v1();
	                                }


                                    if (!found && !segments.empty()) { // remove first segment, and try again!
                                        cout << "Removing a segment while trying to find relation to existing row paths, in parentFace = " << segments[0].parentFace << endl;
                                        segments.erase(segments.begin()); // remove first element
                                    }
                                }
                                if (segments.empty()) 
                                    goto emtpy_segments;


                                // --- find position_i
                                double min_d = 1e7;
                                for (int i = side_finding_i; i <= options; ++i) {
                                    int pi = i;
                                    double d;
                                    if (at_beginning)
                                        d = (segments[0].v1 - current_row_paths.paths[pi].v_first()).squaredNorm();
                                    else
                                        d = (segments[0].v0 - current_row_paths.paths[pi].v_last()).squaredNorm();
                                    if (d < min_d) {
                                        min_d = d;
                                        position_i = i;
                                    }
                                }
                            }
                        }
                    }
                    


                    // --- add segments to paths
                    {
                        if (first_in_row) { // if the paths are empty (i.e. these are the first segments of this row)
                            for (int i = 0; i < segments.size(); ++i)
                                current_row_paths.paths[i].segments.push_back(std::make_shared<Paths::Segment>(segments[i]));
                        }

                        else  // if the paths are not empty, then find where the current segments need to be added (i.e. at the beginning or at the end of each existing path)
                        {
                            bool snap_points_together = true;

                            bool at_the_beginning = segments[0].parentEdges[1] == current_row_paths.paths[0 + position_i].ei_first();
                            bool at_the_end       = segments[0].parentEdges[0] == current_row_paths.paths[0 + position_i].ei_last();
                            if (!at_the_beginning && !at_the_end)
                            {
                                at_the_beginning = are_neighboring_faces(fj, current_row_paths.paths[0 + position_i].segments[0]->parentFace, FE);
                                at_the_end       = are_neighboring_faces(fj, current_row_paths.paths[0 + position_i].segments.back()->parentFace, FE);
                            }
                            if (!at_the_beginning && !at_the_end)
                            {
	                            cerr << "Could not find where to put current segments with fj = " << fj << ". Will skip face." << endl;
                                goto emtpy_segments;
                            }


                            // --- Add segments at the BEGINNING  of the existing paths
                            if (at_the_beginning) 
                            {
                                for (int i = 0; i < segments.size(); ++i) {
                                    int pi = i + position_i;
                                    assert(pi < current_nP);

                                    // check distance
                                    double d = (segments[i].v1 - current_row_paths.paths[pi].v_first()).norm();
                                    if (d > d_threshold) {
                                        cerr << " --- Attention, error! Too large distance between segments: " << d << "for  fj = " << fj << endl; //throw;
                                    }

                                    if (snap_points_together)
                                        segments[i].v1 = current_row_paths.paths[pi].v_first(); // snap the two vertices together

                                    current_row_paths.paths[pi].segments.insert(current_row_paths.paths[pi].segments.begin(), std::make_shared<Paths::Segment>(segments[i]));
                                }
                            }

                            // --- Add segments at the END of the existing paths
                            else 
                            {
                            	for (int i = 0; i < segments.size(); ++i) {
                                    int pi = i + position_i;
                                    assert(pi < current_nP);

                                    // check distance
                                    double d = (segments[i].v0 - current_row_paths.paths[pi].v_last()).norm();
                                    if (d > d_threshold) {
                                        cerr << "--- Attention, error! Too large distance between segments: " << d << " for  fj = " << fj << endl; //throw;
                                    }

                                    if (snap_points_together)
                                        segments[i].v0 = current_row_paths.paths[pi].v_last(); // snap the two vertices together

                                    current_row_paths.paths[pi].segments.push_back(std::make_shared<Paths::Segment>(segments[i]));
                                }
                            }
                        }
                    }

				}

                emtpy_segments:;

                // KEEP COMMENTED OUT, UNLESS DEBUGGING --->
                 //current_row_paths.write_to_file("C:/dev/np3dp_v2/data/costa_minimal_surface_2/output/current_row_paths.txt"); // save current paths
                // KEEP COMMENTED OUT, UNLESS DEBUGGING <---
            }


            // --- add to queue row direct neighbor(s) in the same row
            for (int m=0; m<D[fj]; ++m) {
	            int ei = FE(fj, m);
                if (Emap[ei] == 1) {
	                int nfj = EF(ei, 0) == fj ? EF(ei, 1) : EF(ei, 0);
                    if (nfj >= 0) {
                        if (!visited[nfj]) {
                            R.push_back(nfj);
                        }
                    }
                }
            }
            continue_row:;
        }


        // --------------------------- finished row R
        if (!current_row_paths.check_edge_links())
        {
            cerr << "Wrong edge links in current row!" << endl;
        }
        
        
        // --- append new paths in the correct position, in relation to all existing paths (either at the end, or at the beginning)
        if (paths_collection.paths.empty()) { // if paths are empty, then it doesn't matter where to put them
            paths_collection.insert_paths_at_end(current_row_paths);
        }
        else {
            int nfi, ei;
            bool found_Q = PathTracingHelpers::find_visited_direct_neighbor(fi, 0, visited, Emap, D, FE, EF, nfi, ei);
            if (!found_Q) continue;
            RowVector3d n_dir = guiding_dirs.row(nfi);
            assert(found_Q && 0.9 < n_dir.squaredNorm() < 1.1);

            RowVector3d cen_fi(0, 0, 0), cen_nfi(0, 0, 0); // find the centers of the two faces
            for (int m = 0; m < D[fi]; ++m)  cen_fi += V.row(F(fi, m)) / double(D[fi]);
            for (int m = 0; m < D[nfi]; ++m) cen_nfi += V.row(F(nfi, m)) / double(D[nfi]);

            if ((cen_fi - cen_nfi).dot(n_dir) < 0)
                paths_collection.insert_paths_at_beginning(current_row_paths);
            else
                paths_collection.insert_paths_at_end(current_row_paths);
        }

        
    	// KEEP COMMENTED OUT, UNLESS DEBUGGING --->
         //paths_collection.write_to_file("C:/dev/np3dp_v2/data/costa_minimal_surface_2/output/all_temp_paths.txt"); // save current paths
    	// KEEP COMMENTED OUT, UNLESS DEBUGGING <---
        

        // Get next rows (i.e. direct neighbors of fi that are not on the same row) and add them to Q
        for (int ff : faces_in_current_row)
        {
            for (int m = 0; m < D[ff]; ++m) {
                int ei = FE(ff, m);
                if (Emap[ei] == 0) {
                    int nfi = EF(ei, 0) == ff ? EF(ei, 1) : EF(ei, 0);
                    if (nfi >= 0) {
                        if (!visited[nfi]) {
                            Q.push_back(nfi);
                        }
                    }
                }
            }
        }

    }

    // --- check that all faces where visited
    for (int fi=0; fi<visited.size(); ++fi) {
	    if (visited[fi] == 0) {
            cerr << endl << " --- Attention! Face fi = " << fi << " was not visited!" << endl;
            //throw;
	    }
    }


	// --- remove any paths that have less than 2 segments
	int count = 0;
    for (int pi = paths_collection.paths.size() - 1; pi >= 0 ; --pi) {
        auto& path = paths_collection.paths[pi];
        if (path.segments.size() < 2) {
            paths_collection.paths.erase(paths_collection.paths.begin() + pi);
            ++count;
        }
    }
    cout << "Removed " << count << " paths due to very small size" << endl;


    // --- Label closed paths
    for (int pi=0; pi<paths_collection.paths.size(); ++pi)
    {
        Paths::Path& path = paths_collection.paths[pi];
        if (path.segments[0]->parentEdges[0] == path.segments.back()->parentEdges[1])
            path.is_closed = true;
    }

    // --- if the paths are closed, then align them based on some random subdominant edge sequence
	{
		bool there_are_closed_paths = false;
    	for (const Paths::Path& p : paths_collection.paths){
    		if (p.is_closed){
    			there_are_closed_paths = true;
    			break;
    		}
    	}

    	if (there_are_closed_paths)
    	{
    		// get vertexDegreeQuad and boundaryVerticesQuad
    		VectorXi vertexDegree = QuadMesh::mesh_vertex_degrees(V, F, D);
    		VectorXi boundaryVertices = QuadMesh::vertices_on_boundaries(V, EF, EV);

    		// find a non-boundary face
    		auto count_boundary_vis_in_face = [](int fi, const MatrixXi& F, const VectorXi& boundaryVerticesQuad)->int {
    			int count = 0;
    			for (int m=0; m<F.cols(); ++m){
    				int vi = F(fi, m);
    				if (vi >= 0)
    					if (boundaryVerticesQuad[vi] == 1)
    						++count;
    			}
    			return count;
    		};

    		int fi = 0;
    		while (D[fi] != 4 || count_boundary_vis_in_face(fi, F, boundaryVertices) > 0 || PathTracingHelpers::face_has_rib(fi, D, ERibMap, FE)){
    			++fi;
                if (fi == F.rows()) break;
    		}

    		// get halfedges of a subdominant edge
    		int ei;
    		for (int m = 0; m < D[fi]; ++m) {
    			ei = FE(fi, m);
    			if (Emap[ei] == 1)
    				break;
    		}
    		int he1 = EH(ei,0);
    		int he2 = nextH[twinH[nextH[twinH[he1]]]];

    		// expand cutting from halfedges
    		VectorXi C; vector<int> passed_vis;
    		QuadMesh::expand_cutting_from_he(he1, C, passed_vis, EV.rows(), HE, HV, nextH, twinH, vertexDegree, boundaryVertices, false, true);
    		QuadMesh::expand_cutting_from_he(he2, C, passed_vis, EV.rows(), HE, HV, nextH, twinH, vertexDegree, boundaryVertices, false, true);

    		// update closed paths to start from this edge sequence
            for (Paths::Path& p : paths_collection.paths) {
                if (p.is_closed)
                {
                    int found_si = -1;
	                for (int si=0; si<p.segments.size(); ++si)
	                {
                        int e1 = p.segments[si]->parentEdges[0];
                        // int e2 = p.segments[si]->parentEdges[1];
                        if (C[e1] == 1)
                            found_si = si;
	                }

                    if (found_si >= 0){
	                      std::rotate(p.segments.begin(), p.segments.begin() + found_si, p.segments.end());
                    }
                    else{
                        cout << "!!! DID NOT FIND si FOR UPDATING SEGMENTS ORDER IN CLOSED PATH" << endl;
                    }

                }
            }
    	}
	}
}



void PathsTracing::assign_number_of_paths_per_strip(Piece& p, int default_paths_density, double mesh_scale_factor, double desired_avg_layer_height)
{
	// get all strips
    MatrixXi strips; // Strips x F (with values strips(si, fi) = 1/0 if fi belongs/doesn't belong to strip si)
    VectorXi strips_directions; // Strips x 1 (with values 0/1 for strip in the dominant/subdominant direction)
    vector<array<vector<int>, 2>> strips_vis;
    vector<bool> strips_are_closed;
    vector<vector<int>> StoEis, StoFis; // Strip to Edge indices. For each strip: a vector of ordered edge indices
    MatrixXi FtoS; // Fx2: face to strip (for each strip-direction)
    vector<vector<set<int>>> VtoS; // Vx2: vertex to strip (for each strip-direction):  VtoS[vi][strip_direction][si]

    MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
	hedra::polygonal_edge_topology(p.D, p.F, EV, FE, EF, EFi, FEs, innerEdges);
    VectorXi boundaryFaces = QuadMesh::faces_on_boundaries(p.F, EF);

    Strips::get_all_strips(p.V, p.F, p.D, boundaryFaces, p.Emap, strips, strips_directions,
        strips_vis, strips_are_closed, StoEis, StoFis, FtoS, VtoS, { 0 });

    p.nP_per_strip.setZero(p.F.rows());


    for (int si=0; si<strips.rows(); ++si)
    {
        // get lengths of subdominant edges
        vector<int> eis = StoEis[si];
        vector<double> lengths;

        for (int ei : eis)
	        if (p.Emap[ei] == 1)
                lengths.push_back((p.V.row(EV(ei, 0)) - p.V.row(EV(ei, 1))).norm());

        // get min-max-avg of lengths
        double min = *min_element(lengths.begin(), lengths.end());
        double max = *max_element(lengths.begin(), lengths.end());
        double avg = std::accumulate(lengths.begin(), lengths.end(), 0.0) / lengths.size();
        

        // decide on final number
        double avg_h =  desired_avg_layer_height / mesh_scale_factor;

        int N = std::ceil(0.85 * avg / avg_h);

        if (max / N > 2.5 * avg_h){
            while (max / N > 2.5 * avg_h)
                ++N;
        }
        else if (max / N < 0.5 * avg_h){
        	while (max / N < 0.5 * avg_h)
                --N;
        }


        vector<int> fis = StoFis[si];

        // assign values to the faces of the strip
        for (int fi : fis) {
            p.nP_per_strip[fi] = N;
        }
    }

    // check that all faces have an assigned value
    for (int fi=0; fi<p.F.rows(); ++fi){
	    if (p.nP_per_strip[fi] == 0){
            cout << "--- Attention! fi = " << fi << " does not have a paths-per-strip value. Assigning default." << endl;
            p.nP_per_strip[fi] = default_paths_density;
	    }
    }
}


//////////////////////////////////////
/// --- Single face
/////////////////////////////////////


vector<Paths::Segment> PathsTracing::trace_segments_on_quad_face(int fi, int density, const MatrixXd& V, const MatrixXi& FH, const VectorXi& HE, const VectorXi& nextH,
    const VectorXi& HV, const VectorXi& Emap)
{
    int he = FH(fi, 0);
    int start_he = he;

    while (Emap[HE[he]] != 1 ) { // find a sub-dominant edge in the face
        he = nextH[he];

        if (start_he == he) // if there is no subdominant edge in this face
        { // then take the next of a dominant edge
            while (Emap[HE[he]] != 0)
            {
                he = nextH[he];
                if (start_he == he){cerr << "Could not find any dominant edge in face fi = " << fi << ". Will not trace paths on that face" << endl; return vector<Paths::Segment>();}
            }
            he = nextH[he];
            break;
        }
    }

    vector<int> dominant_eis(2, -1), subdominant_eis(2, -1);
    MatrixXd Vstarts, Vends; // coordinates of edge vertices (start and end)
    Vstarts.setZero(2, 3); Vends.setZero(2, 3);

    Vstarts.row(0) = V.row(HV[he]);
    assert(Emap[HE[he]] == 1 || Emap[HE[he]] == -1); // subdominant or boundary
    subdominant_eis[0] = HE[he];

    he = nextH[he];
    assert(Emap[HE[he]] == 0 || Emap[HE[he]] == -1); // dominant or boundary
    Vends.row(0) = V.row(HV[he]);
    dominant_eis[0] = HE[he];

    he = nextH[he];
    assert(Emap[HE[he]] == 1 || Emap[HE[he]] == -1); // subdominant or boundary
    Vends.row(1) = V.row(HV[he]);
    subdominant_eis[1] = HE[he];

    he = nextH[he];
    assert(Emap[HE[he]] == 0 || Emap[HE[he]] == -1); // dominant or boundary
    Vstarts.row(1) = V.row(HV[he]);
    dominant_eis[1] = HE[he];

    assert(dominant_eis[0] >= 0 && dominant_eis[1] >= 0 && dominant_eis[0] != dominant_eis[1]); // make sure all eis were found
    assert(subdominant_eis[0] >= 0 && subdominant_eis[1] >= 0 && subdominant_eis[0] != subdominant_eis[1]); // make sure all eis were found

    // --- create segments
    vector<Paths::Segment> segments;
    for (int k = 0; k < density; ++k)
    {
        segments.push_back(Paths::Segment());
        Paths::Segment& s = segments.back();

        double p = double(k) / double(density) + 0.5 / double(density);

        s.v0 = p * Vstarts.row(0) + (1. - p) * Vends.row(0);
        s.parentEdges(0) = subdominant_eis[0];

        // End of contour segment: linear combination of v2 and v considering the scalar of v2: t[e2](fi, j)
        s.v1 = p * Vstarts.row(1) + (1. - p) * Vends.row(1);
        s.parentEdges(1) = subdominant_eis[1];

        assert((s.v0 - s.v1).squaredNorm() >= 1e-16); // the segment should never have zero length
        s.parentFace = fi;
    }
    return segments;
}


vector<Paths::Segment> PathsTracing::trace_segments_on_irregular_face(int fi, int density, const MatrixXd& V, const MatrixXi& F, const VectorXi& D, const VectorXi& Emap,
																      const MatrixXd& UV, const MatrixXi& FH, const VectorXi& HE, const VectorXi& prevH, const MatrixXi& EV,
																      const RowVector3d& rot_guiding_dir, double avg_edge_len)
{
    VectorXd spinning_forms(D[fi]); // store the 'spinning form' (i.e. uv change) of each edge
    VectorXd u(D[fi]);

    // count subdominant edges
    int count_subdominant = 0;
    for (int m = 0; m < D[fi]; ++m) {
        int ei = HE[FH(fi, m)];

        if (Emap[ei] == 1)
            ++count_subdominant;
    }

    double e = 0.01; // arbitrary value, might need to be updated..

    VectorXi spinning_form_uncertainty; spinning_form_uncertainty.setZero(D[fi]); // remember how unsure we are for every calculation

    for (int m = 0; m < D[fi]; ++m) {
        double u0 = UV(F(fi, m), 0);
        double u1 = UV(F(fi, (m + 1) % D[fi]), 0);
        int v0 = F(fi, m);
        int v1 = F(fi, (m + 1) % D[fi]);
        RowVector3d edge_vec = V.row(v1) - V.row(v0);
        int sign = rot_guiding_dir.dot(edge_vec) > 0 ? 1 : -1;
        int ei = HE[FH(fi, m)];

        if (Emap[ei] == 0) { // all dominant edges have spinning form 0
            spinning_forms[m] = 0;
        }
    	else if (PathTracingHelpers::is_unit_edge(fi, m, EV, FH, Emap, prevH, HE, D)) { // all unit norm edges have spinning form 1
            spinning_forms[m] = 1 * sign;
    	}
        else 
        { // subdominant or boundary edges can have spinning form 0 <= s <= 1
            double s = abs(abs(u1) - abs(u0)); // only consider absolute values
            if (s < e || s > 1-e) { // all edges should be unit length or less
                // int sign = guiding_dir.dot(edge_vec) > 0 ? 1 : -1; //s / abs(s);
                double decimals = abs(abs(s) - floor(abs(s))); // get rid of additional 'period jumps'
                assert(decimals < 1.0);
                if (decimals < e || decimals > 1-e) // this should result in either 0 or +-1
                {
                    if ((D[fi] >= 4 && count_subdominant == 2 && Emap[ei] == -1) || // if we are on the boundary of a quad with two (complete) subdominant edges
                        edge_vec.norm() < 0.1 * avg_edge_len) // or if the edge length is really small 
						s = 0;
                    else
                        s = 1;
                }
                else
                {
                    s = decimals;
                }
            }
            
            spinning_forms[m] = s * sign;

            // also remember uncertainty
            if (Emap[ei] == -1)
                spinning_form_uncertainty[m] = 2;
            else
                spinning_form_uncertainty[m] = 1;
    	} 
    }

    // --- get starting u (find a corner from the uncut quad mesh, where the UV coords should be integer values)
    int start_m = -1;
    for (int m=0; m<D[fi]; ++m) {
        if (PathTracingHelpers::is_uncut_vertex(fi, m, FH, Emap, prevH, HE)) {
            start_m = m;
            u[m] = 0.0; //round(UV(F(fi, m), 0));
            break;
        }
    }
    if (start_m == -1) { // if start_m was not found
	    start_m = 0;
        u[0] = UV(F(fi, 0), 0);
    }

    double sum = spinning_forms.sum();
    if (abs(sum) > e) {
        cerr << endl << " --- Attention! Sum of spinning forms = " << sum << " for face fi = " << fi << ". Will be updated to equal zero." << endl;
        //cout << "spinning_forms : "; for (int m = 0; m < D[fi]; ++m) cout << spinning_forms[m] << " , "; cout << endl;
        // throw;

        // just update spinning form with highest uncertainty so that the spinning form equals 0
        int max_i; // index with highest uncertainty
        spinning_form_uncertainty.maxCoeff(&max_i);
        double partial_sum = 0;
        for (int m = 0; m < D[fi]; ++m)
            if (m != max_i)
                partial_sum += spinning_forms[m];
        spinning_forms[max_i] = -partial_sum;
        assert(abs(spinning_forms.sum()) < 1e-10);
    }


    // --- get u (i.e. the u coordinate of uv, getting rid of any period jumps)
    for (int m = 1; m < D[fi]; ++m) {
        int i = (start_m + m) % D[fi];
        int i_prev = (start_m + m - 1) % D[fi];
        u[i] = u[i_prev] + spinning_forms[i_prev];

        if (PathTracingHelpers::is_uncut_vertex(fi, i, FH, Emap, prevH, HE))
            u[i] = round(u[i]);
    }




    double u_min = floor(u.minCoeff());
    double iso = u_min + 0.5 / double(density);

    vector<Paths::Segment> segments;
    while (iso < u.maxCoeff())
    {
        Intersections inters;

        // --- find intersections
        for (int m = 0; m < D[fi]; ++m) {
            int vi = F(fi, m);
            int vj = F(fi, (m + 1) % D[fi]);

            double c1 = u[m];
            double c2 = u[(m + 1) % D[fi]];

            double coeff = PathTracingHelpers::is_edge_intersected(c1, c2, iso);

            if (std::isfinite(coeff)) { // if there is an intersection
                int he = FH(fi, m);
                int ei = HE(he);

                Vector3d point = (1 - coeff) * V.row(vi) + coeff * V.row(vj);
                inters.add_intersection(ei, point);
            }
        }

        // --- create segments
        if (inters.is_complete())
        {
            segments.push_back(Paths::Segment());
            Paths::Segment& s = segments.back();

            s.v0 = inters.p1;
            s.parentEdges(0) = inters.e1;
            s.v1 = inters.p2;
            s.parentEdges(1) = inters.e2;
            s.parentFace = fi;
            s.forward_dir = Vector3d(0, 0, 0); // set to zero vector for now
        }

        iso += 1.0 / double(density); // increase isoline value
    }


    //// make sure segments look in the correct direction, with respect to face orientation
    //TODO: attention! This doesn't seem to work! investigate
    if (segments.size() > 1) // Attention! This doesn't work for a face where there's less that 2 segments traced
    {
        const auto& seg0 = segments[0];
        int e = seg0.parentEdges[0];
        Vector3d v0 = seg0.v0;
        Vector3d e_v0 = V.row(EV(e, 0));
        Vector3d e_v1 = V.row(EV(e, 1));
        if ((v0 - e_v1).squaredNorm() < (v0 - e_v0).squaredNorm()) {  // then you need to flip segments direction
            for (auto& seg : segments)
                seg.switch_v0_v1();
        }
    }

    // --- make sure all segments are consistently oriented
    if (segments.size() > 1) {
        for (int si=1; si<segments.size(); ++si) {
            const auto& prev_seg = segments[si-1];
            auto& seg = segments[si];
            if ((prev_seg.v1 - prev_seg.v0).dot(seg.v1 - seg.v0) < 0)
                seg.switch_v0_v1();
        }
    }

    return segments;
}



//////////////////////////////////////
/// --- Helpers
//////////////////////////////////////


void PathsTracing::cleanup(Paths::PathCollection& pathCollection) {
    vector<int> to_remove;
    for (int i = 0; i < pathCollection.paths.size(); ++i)
        if (pathCollection.paths[i].total_length() < 0.02)
            to_remove.push_back(i);

    std::reverse(to_remove.begin(), to_remove.end());
    for (auto& pi : to_remove)
        pathCollection.remove(pi);

    cout << "Cleaning up. Removed " << to_remove.size() << " paths that were too short." << endl;
}

