#include "strips.h"

#include <directional/polygonal_edge_topology.h>

#include "quad_mesh.h"

#define INF 0x3f3f3f3f


void Strips::get_all_strips(const MatrixXd& V, const MatrixXi& F, const VectorXi& D, const VectorXi& boundaryFaces, const VectorXi& Emap,
		MatrixXi& strips, VectorXi& strips_directions, vector<array<vector<int>, 2>>& strips_vis, vector<bool>& strips_are_closed, vector<vector<int>>& eis,
		vector<vector<int>>& fis, MatrixXi& FtoS, vector<vector<set<int>>>& VtoS, const vector<int>& directions)
{
    FtoS.setConstant(F.rows(), 2, -1);
    VtoS.clear(); VtoS.resize(V.rows()); for (int vi = 0; vi < V.rows(); ++vi) { VtoS[vi].resize(2); } // for every vertex, for every strip direction

    // clear any existing values in vectors
    strips.setZero(0, 0);
    strips_directions.setZero(0);
    strips_vis.clear();
    strips_are_closed.clear();
    eis.clear();
    fis.clear();

    // find all strips + strips data in both directions, create maps FtoS and VtoS
    for (int fi = 0; fi < F.rows(); ++fi)
    {
        for (const int strip_direction: directions) {
            if (FtoS(fi, strip_direction) == -1) { // if this strip has not been found
                //cout << endl; for (int i=0; i<D[fi]; ++i) cout << F(fi, i) << endl;

                VectorXi strip; array<vector<int>, 2> strip_vis; bool strip_is_closed; vector<int> strip_eis; vector<int> strip_fis;
                bool success = Strips::get_strip_from_fi(D, F, boundaryFaces, Emap, fi, strip_direction, strip, strip_vis, strip_is_closed, strip_eis, strip_fis);

                // --- check result
                if (success){ // make sure that fi is in the strip
                    if (strip[fi] != 1){
	                    cerr << "While creating all strips, found a strip that does not contain the original fi = " << fi << endl;
                        success = false;
                    }
                }
                if (success) { // make sure that none of the faces of the strip have been used
                    for (int fj = 0; fj < F.rows(); ++fj) {
                        if (strip[fj]) {
                            if (FtoS(fj, strip_direction) != -1) {
                                cerr << "While creating all strips, found a new strip that uses an already occupied face fj = " << fj << endl;
                                success = false;
                                break;
                            }
                        }
                    }
                }

                // --- save strip
                if (success) {
                    int si = strips.rows();
                    strips.conservativeResize(si + 1, F.rows());
                    strips.row(si) = strip.transpose();

                    strips_directions.conservativeResize(si + 1);
                    strips_directions[si] = strip_direction;

                    for (int fi = 0; fi < F.rows(); ++fi) // for each face that belongs to the strip, s
                        if (strip[fi]){
                            FtoS(fi, strip_direction) = si;
                            //cout << fi << endl;
                        }

                    for (int k = 0; k < strip_vis.size(); ++k) {
                        for (int vi : strip_vis[k]) {
                            VtoS[vi][strip_direction].insert(si);
                        }
                    }

                    strips_vis.push_back(strip_vis);
                    strips_are_closed.push_back(strip_is_closed);
                    eis.push_back(strip_eis);
                    fis.push_back(strip_fis);

                    //Helpers::write_matrix_to_txt(strips, "C:/dev/np3dp_v2/data/square_gyroid/output/temp_strips.txt");
                    //Helpers::write_vector_to_txt(strips_directions, "C:/dev/np3dp_v2/data/square_gyroid/output/temp_strips_directions.txt");
                    //cout << endl << endl;
                }
            }
        }
    }
}



namespace GraphHelpers
{
    GraphVertex merge_vertices(const Np3dpContext* np3dp, StripGraph& G, vector<GraphVertex>& vis, int common_vi) // returns the vertex descriptor of the new vertex (that contains all the merged vertices)
    { // common_vi = -1 when no common_vi exists
        // --- create new merged vertex
	    GraphVertex g_merged = boost::add_vertex(G);
        G[g_merged].color = "black";
        G[g_merged].shape = "invtriangle";

        set<Edge> edges_to_remove;

	    for (GraphVertex g : vis){
            if (g != common_vi) // if it's the common vi, then ignore 
            { 
                for (GraphVertex neighbor : boost::make_iterator_range(boost::adjacent_vertices(g, G))) { // iterate through all edges of g
                    Edge ei_to_remove = boost::edge(g, neighbor, G).first;
                    edges_to_remove.insert(ei_to_remove);

                    if (!Helpers::in_vector(vis, neighbor)) { // add_new_edge between g_merged and the other neighbor (if neighbor is not going to be merged)
                        // get replacing_si for new edge
                        int r_si;
                        if (get(&GraphEdgeProps::replacing_strip, G)[ei_to_remove] == -1) // if ei_to_remove is not an already merged edge
                            r_si = G[g].si;
                        else
                            r_si = get(&GraphEdgeProps::replacing_strip, G)[ei_to_remove];

                        // check if edge already exists with the same r_si
                        bool identical_edge_exists = false;
                        auto e_pair = boost::edge(g_merged, neighbor, G);
                        if (e_pair.second)
                            if (get(&GraphEdgeProps::replacing_strip, G)[e_pair.first] == r_si)
                                identical_edge_exists = true;

                        if (!identical_edge_exists)
                        {
                            Edge ei = boost::add_edge(g_merged, neighbor, G).first; // add edge towards merged vertex
                            get(&GraphEdgeProps::replacing_strip, G)[ei] = r_si;
                            //StripsTopology::write_graphviz(np3dp, DATA_PATH + np3dp->output_folder + "merged_sings_cD.dot", DATA_PATH + np3dp->output_folder + "_merged_sings_cD.png", G);
                        }
                    }
                }
            }
	    }


        // --- remove edges
        for (std::set<Edge>::iterator it = edges_to_remove.begin(); it != edges_to_remove.end(); ++it) {
            const Edge ei = *it;
            boost::remove_edge(ei, G); // delete existing edge
        }
        //StripsTopology::write_graphviz(np3dp, DATA_PATH + np3dp->output_folder + "merged_sings_cD.dot", DATA_PATH + np3dp->output_folder + "_merged_sings_cD.png", G);


        // --- if there is a common_vi, then add edge between g_merged and common_vi
        if (common_vi >= 0)
            boost::add_edge(g_merged, common_vi, G);


        // ---  remove vertices
        sort(vis.begin(), vis.end(), [](GraphVertex a, GraphVertex b)->bool {return a > b;}); // sort from highest to lowest index
    	for (const GraphVertex vi : vis) {
            if (vi != common_vi) // if it's the common vi, then ignore 
            {
                boost::remove_vertex(vi, G); // delete existing edge
                if (vi < g_merged)
                    --g_merged; // update vertex index
            }
        }

        return g_merged;
    }

    void convert_strip_indices_to_graph_indices(const Np3dpContext* np3dp, const StripGraph& G, const set<int>& strips, vector<GraphVertex>& graph_vertices)
    {
    	const QuadMesh& q = np3dp->quad_mesh;
        for (std::set<int>::iterator it = strips.begin(); it != strips.end(); ++it) {
            const int si = *it;
            GraphVertex gi = q.strip_networks.StoG[si];
            graph_vertices.push_back(gi);
            if (G[gi].si != si) throw invalid_argument("Wrong selection of graph or of graphVertex. cD[gi].si != si"); // make sure that the correct graph vertices were found
        }
    }

    bool two_strips_are_neighboring(const Np3dpContext* np3dp, int s1,int s2)
    {
        const QuadMesh& q = np3dp->quad_mesh;
	    const std::array<vector<int>, 2>& strip_vis_1 = q.strip_networks.strip_index_to_properties(s1).strip_vis;
        const std::array<vector<int>, 2>& strip_vis_2 = q.strip_networks.strip_index_to_properties(s2).strip_vis;

        auto common_elements_exist = [](const vector<int>& vis1, const vector<int>& vis2)->bool
        {
	        for (int v1 : vis1)
		        for (int v2 : vis2)
                    if (v1 == v2)
                        return true;
	        return false;
        };

        for (int side=0; side<2; ++side){
            const vector<int>& vis1 = strip_vis_1[side];
            for (int other_side = 0; other_side <2; ++other_side){
                const vector<int>& vis2 = strip_vis_2[other_side];
                if (common_elements_exist(vis1, vis2))
                    return true;
            }
        }
        return false;
    }

    vector<Edge> get_vertex_edges(const StripGraph& G, int vi)
    {
        vector<Edge> edges;
        typename GraphTraits::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(vi, G); ei != ei_end; ++ei) {
            //auto source = boost::source(*ei, G);
            //auto target = boost::target(*ei, G);
            //std::cout << "There is an edge from " << source << " to " << target << std::endl;
            edges.push_back(*ei);
        }
        return edges;
    }
}


void Strips::create_strip_graphs(Np3dpContext* np3dp)
{
    QuadMesh& q = np3dp->quad_mesh;
    StripNetworks& sn = q.strip_networks;
    sn.invalidate();

    StripGraph& D = sn.Dgraph;
    StripGraph& cD = sn.cDgraph;
    MatrixXi& FtoS = sn.FtoS;
    VectorXi& StoG = sn.StoG;
    VectorXi& StoD = sn.StoD;
    VectorXi& EtoS = sn.EtoS;
    MatrixXi& strips = sn.strips;
    vector<vector<int>>& StoEis = sn.StoEis;
    vector<vector<int>>& StoFis = sn.StoFis;
    vector<vector<set<int>>>& VtoS = sn.VtoS;

    VectorXi strips_directions; // Strips x 1 (with values 0/1 for strip in the dominant/subdominant direction)
    vector<array<vector<int>, 2>> strips_vis;
    vector<bool> strips_are_closed;
    get_all_strips(q.V, q.F, q.D, q.boundaryFaces, q.Emap, strips, strips_directions, strips_vis, strips_are_closed, StoEis, StoFis, FtoS, VtoS);
    
    Helpers::write_matrix_to_txt(strips, DATA_PATH + np3dp->output_folder + "strips.txt");
    Helpers::write_vector_to_txt(strips_directions, DATA_PATH + np3dp->output_folder + "strips_directions.txt");

    int nS = strips.rows(); // nS is the sum of strips in dir 0 and dir 1

    StoG.setConstant(nS, 1, -1);
    StoD.setConstant(nS, 1, -1);

    // --- create graph nodes. Each node corresponds to one strip.
    for (int si = 0; si < nS; ++si) {
        int strip_dir = strips_directions[si];
        StripGraph& G = strip_dir == 0 ? D : cD;

        auto v = boost::add_vertex(G);
        G[v].strip = strips.row(si);
        G[v].si = si;
        G[v].strip_is_closed = strips_are_closed[si];
        G[v].strip_vis = strips_vis[si];

        StoG[si] = v;
        StoD[si] = strip_dir;
    }

    // --- create graph edges
    for (int k = 0; k < 2; ++k)
    {
        StripGraph& G = k == 0 ? D : cD;
        const int strip_dir = k == 0 ? 0 : 1;

        for (int vi = 0; vi < boost::num_vertices(G); ++vi)
        {
            //const VectorXi& strip = get(boost::vertex_attribute, G)[vi].strip;
            const VectorXi& strip = G[vi].strip;

            for (int fi = 0; fi < q.F.rows(); ++fi) { // go through all the faces that belong to the strip
                if (strip[fi]) {
                    for (int j = 0; j < q.D[fi]; ++j) { // find fi neighbors that do not belong in the same strip
                        int ei = q.FE(fi, j);
                        int nf = q.EF(ei, 0) == fi ? q.EF(ei, 1) : q.EF(ei, 0);
                        if (nf >= 0) {
                            if (!strip[nf]) {
                                int other_si = FtoS(nf, strip_dir); // find which strip the neighbors belong to

                                if (other_si >= 0)
                                {
                                    int other_vi = StoG(other_si);

                                    { // debug check
                                        //const VectorXi& other_strip = get(boost::vertex_attribute, G)[other_vi].strip;
                                        const VectorXi& other_strip = G[other_vi].strip;
                                        const VectorXi& other_strip_from_m = strips.row(other_si);
                                        if ((other_strip - other_strip_from_m).cwiseAbs().sum() != 0) throw invalid_argument("Found the wrong strip!");
                                    }

                                    std::pair<Edge, bool> r = boost::edge(vi, other_vi, G); // check if edge already exists
                                    if (!r.second)
                                        boost::add_edge(vi, other_vi, G); // add graph edge between strip + other strip
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "Created D graph with " << boost::num_vertices(D) << " vertices, and " << boost::num_edges(D) << " edges." << endl;
    std::cout << "Created cD graph with " << boost::num_vertices(cD) << " vertices, and " << boost::num_edges(cD) << " edges." << endl;

    // change cDgraph color to green
    for (int vi = 0; vi < boost::num_vertices(cD); ++vi)
        cD[vi].color = "green";

    // --- get EtoS (edge to strip)
    EtoS.setConstant(q.EV.rows(), -1);
    for (int si = 0; si < StoD.rows(); ++si) {
        for (int ei : StoEis[si]) {
            int dir = StoD[si];
            if (q.Emap[ei] == dir) {
                int v0 = q.EV(ei, 0); int v1 = q.EV(ei, 1);
                cerr << "The edges contained in a strip should have Emap = !bool(dir). Edge vertices : v0 = " << v0 << " , v1 = " << v1 << endl; throw;
            }
            if (EtoS[ei] != -1 && EtoS[ei] != si) {
                int other_si = EtoS[ei];
                int v0 = q.EV(ei, 0); int v1 = q.EV(ei, 1);
                std::cout << "Edge already found in other strip!. Edge vertices : v0 = " << v0 << " , v1 = " << v1 << ". Current si = " << si << ", other si = " << other_si << endl; throw;
            }
            EtoS[ei] = si;
        }
    }

    sn.blockedStrips.setZero(nS); // set all strips as not-blocked

    sn.has_strip_networks_data = true;
}


void Strips::find_strips_sequences_between_two_vertices(Np3dpContext* np3dp, int vi, int vj,
    vector<vector<int>>& counterstrip_sequences, vector<vector<int>>& vis_sequences, vector<bool>& aligned, vector<vector<int>>& strip_sequences, vector<vector<int>>& alignment_mesh_paths)
{
    // finds all sequences of strips in the dominant and subdominant directions. For all counter-strip sequences check alignment
    counterstrip_sequences.clear();
    strip_sequences.clear();
    alignment_mesh_paths.clear();

    const VectorXi& Emap = np3dp->quad_mesh.Emap;


    // --- Get all sequences of counter-strips
    Strips::get_all_counterstrip_sequences_between_two_vertices(np3dp, vi, vj, 0, Emap, counterstrip_sequences, vis_sequences);


    // --- Check if there are alignments between two vertices along all the found sequences
    int direction = 0; // always consider alignment in the dominant direction
    int other_direction = 1;
    aligned = Strips::two_vertices_are_aligned(np3dp, vi, vj, direction, Emap, counterstrip_sequences, vis_sequences, alignment_mesh_paths);

    if (std::all_of(aligned.begin(), aligned.end(), [](bool item)->bool {return item == true; })) {
        cout << "\nThe selected vertices are aligned in all dominant directions!" << endl; // No need to search for intermediate blue stips if they are already aligned in every direction
        strip_sequences.resize(counterstrip_sequences.size());
    }
    else
    {
        // --- create VE lists
        vector<vector<int>> VE(np3dp->quad_mesh.V.rows());
        for (int ei = 0; ei < np3dp->quad_mesh.EV.rows(); ++ei){
	        VE[np3dp->quad_mesh.EV(ei, 0)].push_back(ei);
            VE[np3dp->quad_mesh.EV(ei, 1)].push_back(ei);
        }

        for (int i = 0; i < counterstrip_sequences.size(); ++i)
        {
            vector<int> strip_sequence;
            if (!aligned[i])
            {
                int C_vi = vis_sequences[i].back();


                vector<int> path;
                vector<vector<int>> VV = np3dp->quad_mesh.adjacency_list();
                bool found = np3dp->quad_mesh.shortest_path(np3dp->quad_mesh.selection.selected_vis[1], {C_vi}, VV, path);
                if (!found) {cerr << "Could not find mesh path between second singularity and third corner of the square" << endl; goto ctn; }

                for (int k=0; k< path.size()-1; ++k)
                {
	                int v0 = path[k];
                    int v1 = path[k+1];
                    vector<int> common_elements = Helpers::common_elements_in_vectors(VE[v0], VE[v1]);
                    if (!common_elements.size() == 1){cerr << "No common elements in the VE list for vis :" << v0 << " , " << v1 << endl; throw; }
                    int ei = common_elements[0];
                    int si = np3dp->quad_mesh.strip_networks.EtoS[ei]; // strip index
                    if (np3dp->quad_mesh.strip_networks.StoD[si] == direction)
						strip_sequence.push_back(si);
                }

            	//deque<GraphVertex> path = StripsTopology::shortest_strips_path_between_two_mesh_vertices(np3dp, C_vi, np3dp->selected_vis[1], direction, side_id);
                //for (int i = 0; i < path.size(); ++i) {// remember: path holds graph vertices
                //    const GraphVertex& v = path[i];
                //    strip_sequence.push_back(np3dp->Dgraph[v].si);
                //}
            }
            ctn:;
        	strip_sequences.push_back(strip_sequence);
        }

        if (np3dp->quad_mesh.strip_networks.selected_sequence_index >= counterstrip_sequences.size())
            np3dp->quad_mesh.strip_networks.selected_sequence_index = 0;

        if (counterstrip_sequences.empty()) {
            std::cout << "Could not find any connecting sequences!" << endl;
            return;
        }

        std::cout << "\nCurrent counterstrips sequence_index = 0 " << ", out of " << counterstrip_sequences.size() - 1 << " total indices. Press SPACE to change selected sequence index" << endl;
        std::cout << "The selected vertices are **" << strip_sequences[np3dp->quad_mesh.strip_networks.selected_sequence_index].size() << "** strips apart in the blue direction." << endl;
    }
}


void Strips::find_tracked_alignments(Np3dpContext* np3dp, vector<AlignedVertices>& aligned_vis)
{
    cout << "Checking for tracked vis alignments... ";
    const VectorXi& Emap = np3dp->quad_mesh.Emap;

    int direction = 0; // Only consider alignment in the dominant direction

    // --- find all pairs of singularities that are aligned
    for (int i = 0; i < np3dp->quad_mesh.vis_to_check_for_alignments.size(); ++i)
    {
            int vi = np3dp->quad_mesh.vis_to_check_for_alignments[i].first;
            int vj = np3dp->quad_mesh.vis_to_check_for_alignments[i].second;

            vector<vector<int>> counterstrip_sequences; vector<vector<int>> vis_sequences;
            get_all_counterstrip_sequences_between_two_vertices(np3dp, vi, vj, direction, Emap, counterstrip_sequences, vis_sequences);

        	vector<vector<int>> all_mesh_paths;
            vector<bool> aligned = two_vertices_are_aligned(np3dp, vi, vj, direction, Emap, counterstrip_sequences, vis_sequences, all_mesh_paths);

            if (std::any_of(aligned.begin(), aligned.end(), [](bool item)->bool {return item == true; })) // if any of the routes are aligned
                aligned_vis.emplace_back(vi, vj, direction, aligned, all_mesh_paths, counterstrip_sequences);
    }
     cout << "Checked " << np3dp->quad_mesh.vis_to_check_for_alignments.size() << " pairs, and found " << aligned_vis.size() << " alignments." << endl;
}


int Strips::get_all_counterstrip_sequences_between_two_vertices(const Np3dpContext* np3dp, int vi, int vj, int direction, const VectorXi& Emap, vector<vector<int>>& counterstrip_sequences, vector<vector<int>>& vis_sequences)
{
    // finds all the strips that are in the other_direction
	const int other_direction = static_cast<int>(!static_cast<bool>(direction));


    // --- get counter-strips of two vertices
    // const set<int>& counterstrips_i = np3dp->VtoS[vi][other_direction];
    // const set<int>& counterstrips_j = np3dp->VtoS[vj][other_direction];

    const QuadMesh& q = np3dp->quad_mesh;

	// --- get halfedge datastructure
    auto get_all_outgoing_halfedges_in_direction = [=](int vi, int direction)->vector<int>
    {
        // --- get all outgoing half-edges from vi with the color in *direction*
        vector<int> outgoing_heis;
        {
            int start_he = q.VH[vi];
            int he = start_he;
            do
            {
                if (q.HV[he] != vi) throw invalid_argument("Wrong halfedge");
                if (Emap[q.HE[he]] == direction)
                    outgoing_heis.push_back(he);

                he = q.twinH[he]; if (he == -1) break;
                he = q.nextH[he];
            } while (he != start_he && he != -1);
        }
        return outgoing_heis;
    };

    // --- get all outgoing blue half-edges from vi with the color in *direction*
    vector<int> green_halfedges_vj = get_all_outgoing_halfedges_in_direction(vj, other_direction);
    vector<int> blue_halfedges_vi = get_all_outgoing_halfedges_in_direction(vi, direction);

	counterstrip_sequences.clear();
    counterstrip_sequences.resize(blue_halfedges_vi.size());
    vis_sequences.clear();
    vis_sequences.resize(blue_halfedges_vi.size());

    // --- expand green edges from vj
    VectorXi C;
    for (int he : green_halfedges_vj)
		QuadMesh::expand_cutting_from_he(he, C, vector<int>(), q.EH.rows(), q.HE, q.HV, q.nextH, q.twinH, q.vertexDegree, q.boundaryVertices, false, true);

    // --- walk from each outgoing halfedge until you reach: (a) a boundary, (b) a singularity, (c) an edge with C[ei]=1
    vector<bool> to_remove(blue_halfedges_vi.size(), true);

	for (int i=0; i<blue_halfedges_vi.size(); ++i)
    {
        vis_sequences[i].push_back(vi);
        int he = blue_halfedges_vi[i];
	    vector<int>& sequence = counterstrip_sequences[i];
        int prev_csi = -1;

		VectorXi passed_strips; passed_strips.setZero(q.strip_networks.StoD.rows());
		VectorXi passed_he; passed_he.setZero(q.HE.rows());

        vector<int> passed_vis;

	    while (he >= 0)
	    {
            passed_vis.push_back(q.HV[he]);

            // do not allow halfedges to be passed twice
            if (passed_he[he] == 1) break;
            passed_he[he] = 1;

            // find counterstrip that contains he
            int csi = q.strip_networks.EtoS[q.HE[he]]; // cout << csi << endl;
            if (csi == -1) break; // reached a boundary

            // allow crossing counter-strips up to two times
            if (passed_strips[csi] > 2) break; 
            if (prev_csi != csi){ // only count as an additional pass of the strip, if we have exited the strip and now re-entering it. 
            	passed_strips[csi] += 1;
                prev_csi = csi;

                sequence.push_back(csi);
            }

            int v = q.HV[q.twinH[he]]; // third vertex of current square
            vis_sequences[i].push_back(v);

            /*// if reached counter-strip contained by vj, exit
            if (Helpers::in_set(counterstrips_j, csi)) { 
                if (!(vi == vj && sequence.size() == 1)) {// do not exit if we are trying to find routes from a vertex to itself, and it's the first strip visited
	                to_remove[i] = false;
            		break;
                }
            }*/

            // if reached edge with C[ei] = 1, exit
            int e1 = q.HE[q.nextH[he]];
            if (C[e1] == 1) {to_remove[i] = false; break;}
            int e2; int other_he = q.prevH[q.twinH[he]]; if (other_he >= 0) e2 = q.HE[other_he];
            if (C[e2] == 1) {to_remove[i] = false; break;}


            // if reached a singularity, exit
            if (q.vertexDegree[q.HV[he]] != 4 && sequence.size() != 1){
                to_remove[i] = false; // keep this path? I'm not sure...
	            break;
            }

            // find next he
            he = q.twinH[q.nextH[he]]; if (he == -1) break;
            he = q.nextH[he]; if (he == -1) break;

            if (Emap[q.HE[he]] != direction) break;
	    }

        // cout << endl; for (int vi : passed_vis) cout << vi << endl; cout << endl;
    }

    // --- also remove any identical sequences
    for (int i = 0; i < counterstrip_sequences.size(); ++i)
        for (int j = i + 1; j < counterstrip_sequences.size(); ++j)
            if (counterstrip_sequences[i].size() == counterstrip_sequences[j].size())
                if (Helpers::common_elements_in_vectors(counterstrip_sequences[i], counterstrip_sequences[j]).size() == counterstrip_sequences[i].size())
                    to_remove[i] = true;


    // --- removed flagged sequences
    for (int i= counterstrip_sequences.size()-1; i>=0; --i)
		if (to_remove[i])
		{
			counterstrip_sequences.erase(counterstrip_sequences.begin() + i);
            vis_sequences.erase(vis_sequences.begin() + i);
		}

    return counterstrip_sequences.size();
}


vector<bool> Strips::two_vertices_are_aligned(const Np3dpContext* np3dp, int vi, int vj, int direction, const VectorXi& Emap, const vector<vector<int>>& counterstrip_sequences, const vector<vector<int>>& vis_sequences, vector<vector<int>>& mesh_paths)
{
    const QuadMesh& q = np3dp->quad_mesh;

    // --- measure common strips
    const set<int>& strips_i = q.strip_networks.VtoS[vi][direction];
    const set<int>& strips_j = q.strip_networks.VtoS[vj][direction];
    set<int> intersection;
    set_intersection(strips_i.begin(), strips_i.end(), strips_j.begin(), strips_j.end(), std::inserter(intersection, intersection.begin()));

    int N = counterstrip_sequences.size();
    if (vis_sequences.size() != N) {cerr << "While searching if two vertices are aligned (vi = " << vi << " , vj  " << vj << "), ran into error: C_vis.size() != N" << endl; throw;}


    // --- check if they are aligned across each sequence
    vector<bool> aligned;
	for (const vector<int>& vis_seq : vis_sequences)
		aligned.push_back(vis_seq.back() == vj);


    mesh_paths.resize(N);
    for (int i = 0; i < N; ++i)
    {
        if (aligned[i])
        {
            const vector<int>& strips_path = counterstrip_sequences[i];

            // create mesh adjacency lists only using mesh edges whose vertices belong to the counter-strips of the current path
            vector<vector<int>> VV = q.adjacency_list_with_Emap(direction);
            //Helpers::write_VV_to_txt(DATA_PATH + np3dp->output_folder + "VV.txt", VV);

            // find shortest path between mesh vertices
            bool mesh_path_exists = q.shortest_path(vi, { vj }, VV, mesh_paths[i]);
            if (!mesh_path_exists) {cerr << "Could not find mesh path for vertices that are aligned: vi = " << vi << " , vj = " << vj << endl; throw; }
        }
    }

    return aligned;
}


bool Strips::get_strip_from_fi(const VectorXi& D, const MatrixXi& F , const VectorXi& boundaryFaces, const VectorXi& Emap, int fi, 
		                   int strip_direction, VectorXi& strip, array<vector<int>, 2>& strip_vis, bool& strip_is_closed, vector<int>& strip_eis, 
						   vector<int>& strip_fis, bool update_f0_if_boundary)
{
	// --- get halfedge datastructure
	MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
	hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
    if (EV.rows()!= Emap.rows()) throw invalid_argument("While searching for strip from fi, EV.rows() != Emap.rows()");
    MatrixXi EH, FH; VectorXi VH, HE, HF, HV, nextH, prevH, twinH;
    hedra::dcel(D, F, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

    if (EV.rows()!= Emap.rows()) throw invalid_argument("While searching for strip from fi, EV.rows() != Emap.rows()");

    strip_vis[0].clear(); strip_vis[1].clear();
    strip_is_closed = false;
    strip_eis.clear();
    strip_fis.clear();

    strip.setZero(F.rows());

    const int other_direction = static_cast<int>(!static_cast<bool>(strip_direction));

    // --- get starting face fi
    int original_fi = fi;
    if (boundaryFaces[fi] && update_f0_if_boundary) { // if fi is a boundary face in the direction of the strip, then use a neighboring face as a start
        // if updating original_fi there is a bug; there is a risk that the selected f0 will NOT be part of the final strip
    	if (D[fi] != 4 || QuadMesh::count_Emap_values_of_face_edges(fi, D, FE, Emap, other_direction) <2) {
            // cout << "Starting face of strip (fi = " << fi << ") is on boundary, trying to find a better starting face" << endl;
        	std::deque<int> fis;
        	fis.push_back(fi);

            VectorXi passed; passed.setZero(F.rows()); passed[fi] = 1;
            bool found_non_boundary_face = false;
            while (!fis.empty() && !found_non_boundary_face)
            {
                fi = fis[0]; fis.pop_front();

    			for (int i = 0; i < D[fi]; ++i) {
			        int he = FH(fi, i);
			        int ei = HE[he];
			        if (Emap[ei] == other_direction)
			        {
		                he = twinH[he]; if (he == -1) continue;
		                int nfi = HF[he];

                        if (!passed[nfi]){
                            //cout << endl; for (int i = 0; i < D[nfi]; ++i) cout << Fquad(nfi, i) << endl;

	                        if (!boundaryFaces[nfi]) {
	                            found_non_boundary_face = true;
	                            fi = nfi;
	                            break;
	                        }
	                        fis.push_back(nfi);
                            passed[nfi] = 1;
                        }
			        }
    			}
            }
            if (!found_non_boundary_face)
                fi = original_fi;
            // else
            //     cout << "Original strip fi = " << original_fi << " is on boundary, so we use as starting fi instead the face = " << fi << endl;
        }
    }

    strip[fi] = 1;
    strip_fis.push_back(fi);

    // --- get starting halfedges
    vector<int> starting_halfedges;
    for (int i = 0; i < D[fi]; ++i) {
        int he = FH(fi, i);
        int ei = HE[he];
        if (Emap[ei] == other_direction)
            starting_halfedges.push_back(he);
    }
    if (starting_halfedges.size() != 2){
	    //cout << "Could not find strip" << endl;
    	return false;
    }

    auto add_to_vector = [](int i, vector<int>&v, int item)
	{
        if(!Helpers::in_vector(v, item)){
            if (i == 0) // add at the end
                    v.push_back(item);
            else  // add at the beginning 
	            v.insert(v.begin(), item);
        }
    };

    // --- expand from each starting halfedge
    for (int i = 0; i < starting_halfedges.size(); ++i) {
        int he = starting_halfedges[i];

        vector<int>& vis_A = i == 0 ? strip_vis[0] : strip_vis[1];
        vector<int>& vis_B = i == 0 ? strip_vis[1] : strip_vis[0];

        while (he >= 0)
        {
            add_to_vector(i, strip_eis, HE[he]);
            add_to_vector(i, vis_A, HV[he]);
            
            he = twinH[he]; if (he == -1) break;

            add_to_vector(i, vis_B, HV[he]);

            const int nfi = HF[he];
            if (strip[nfi]){ strip_is_closed = true; break;}
            strip[nfi] = 1; // mark face as part of the strip

            add_to_vector(i, strip_fis, nfi);

            if (D[nfi] == 4) // for regular faces : just go forward two he steps 
            {
                he = nextH[nextH[he]];

                if (twinH[he] == -1) {// if we are on a boundary
                    add_to_vector(i, vis_B, HV[prevH[prevH[prevH[he]]]]);
                }
            }
            else // for irregular faces : go forward until you find the next halfedge with the correct color
            {
                if (D[nfi] == 3) {
                    for (int kk=0; kk<3; ++kk){ // find vi that has not been added to vis_A or vis_B
                        int vi = F(nfi, kk);
	                    if (!Helpers::in_vector(vis_A, vi) && !Helpers::in_vector(vis_B, vi)){
		                    if (vis_A.size() < vis_B.size())
                                add_to_vector(i, vis_A, vi);
                            else
                                add_to_vector(i, vis_B, vi);
                            break;
	                    }
                    }
                    //cout << "vis_A : " << endl; Helpers::print_vector(vis_A);
                    //cout << "vis_B : " << endl; Helpers::print_vector(vis_B);
                    break; // finished strip
                }

                if (D[nfi] >= 5)
                {
                    const int starting_he = he; // save starting value
                    bool found_extra_vi = false;

                    do {
                        he = nextH[he];
                        if (Emap[HE[he]] == other_direction)
                            break;
                        if (!Helpers::in_vector(vis_A, HV[he])) {
                            add_to_vector(i, vis_A, HV[he]);
                            found_extra_vi = true;
                        }

                    } while (he != starting_he);

                    if (!found_extra_vi)
                    {
                        he = starting_he;
                        do {
                            he = prevH[he];
                            if (Emap[HE[he]] == other_direction)
                                break;
                            if (!Helpers::in_vector(vis_B, HV[he])) {
                                add_to_vector(i, vis_B, HV[he]);
                                found_extra_vi = true;
                            }
                        } while (he != starting_he);
                    }

                    if (!found_extra_vi){
	                    cerr << "Could not find extra vi in a D5 face!" << endl;
                    }

                    if (he == starting_he) break;
                }
            }
        }
    }

    // --- consolidate strip_vis into two sorted vectors
	// reverse vis coming from the second halfedge

    /*cout << endl;
	Helpers::print_vector(strip_vis[0]);
	cout << endl;
	Helpers::print_vector(strip_vis[1]);
	cout << endl;*/
    /* //make sure there's no duplicate vertices (Attention! There can be duplicate vertices when the strip wraps around itself)
    for (int vi : strip_vis[0]) 
        if (Helpers::in_vector(strip_vis[1], vi)){
	        cerr << "vi = " << vi << " is both in strip_vis[0] and strip_vis[1]. This should never happen" << endl;
            cout << "strip_vis[0] : " << endl; Helpers::print_vector(strip_vis[0]);
            cout << "strip_vis[1] : " << endl; Helpers::print_vector(strip_vis[1]);
        	throw;
        }*/

    return true;
}


deque<GraphVertex> Strips::shortest_strips_path_between_two_mesh_vertices(const Np3dpContext* np3dp, int v1, int v2, int direction, int side_id)
{
    const QuadMesh& q = np3dp->quad_mesh;
    vector<ShortestPathsResults> results;

    const StripGraph& G = direction == 0 ? q.strip_networks.Dgraph : q.strip_networks.cDgraph;

    // find start and end strips
    const set<int>& strips_start = q.strip_networks.VtoS[v1][direction];
    const set<int>& strips_end = q.strip_networks.VtoS[v2][direction];

	//if (TODO!)
	//{
	//    std::cout << "The two vertices are already aligned in that direction" << endl;
	//    return deque<GraphVertex>(); // return empty path
	//}

    // get adjacency lists in blue or green direction
    VectorXi Emap = q.Emap;
    if (side_id == 1) MeshHelpers::flip_Emap(Emap);
    vector<vector<int>> VV = q.adjacency_list_with_Emap(direction); // only considering edges of the direction

    //std::cout << "Searching shortest path in direction = " << direction << endl;

    for (int si_start : strips_start) {
        GraphVertex start_node = q.strip_networks.StoG[si_start];

        for (int si_end : strips_end) {
            GraphVertex end_node = q.strip_networks.StoG[si_end];

            double dist = 0;
            std::deque<GraphVertex> path;
        	vector<int> distmap, predmap;

            if (Graph::shortest_path(G, start_node, { end_node }, distmap, predmap)) {
                dist = distmap[end_node] + 1;
                //cout << "distance = " << dist << endl;
                
                for (GraphVertex current = end_node; current != G.null_vertex() && predmap[current] != current && current != start_node;) {
                    path.push_front(predmap[current]);
                    current = predmap[current];
                }
                //std::cout << "Path graph nodes from #" << start_node << " to #" << end_node << ": ";
                //std::copy(path.begin(), path.end(), std::ostream_iterator<Vertex>(std::cout, ", "));
                //std::cout << end_node << "\n";
                //std::cout << "Path from strip #" << G[start_node].si << " to #" << G[end_node].si << ": ";
                //for (auto& p : path) std::cout << G[p].si << " , ";
                //std::cout << G[end_node].si << endl;

                path.push_back(end_node);
			}
            else
            {
                cerr << "Could not find shortest paths between strip " << si_start << " and " << si_end << endl; throw;
            }
            results.emplace_back(start_node, end_node, si_start, si_end, dist, path);
        }
    }

    std::sort(results.begin(), results.end(), [](const ShortestPathsResults& a, const ShortestPathsResults& b)->bool {return a.dist < b.dist; });
    //std::cout << "Found shortest path with dist = " << results[0].dist << endl; // " and path : ";
    // std::copy(results[0].path.begin(), results[0].path.end(), std::ostream_iterator<Vertex>(std::cout, ", ")); std::cout << endl;

    return results[0].path;
}


void Strips::write_graphviz(const Np3dpContext* np3dp, const std::string& dot_filename, const std::string& png_filename, StripGraph& G)
{
    const QuadMesh& q = np3dp->quad_mesh;

    // update vertex display information
    for (int vi=0; vi<boost::num_vertices(G); ++vi) {
        int si = G[vi].si;
	    G[vi].label = std::to_string(si);
        G[vi].label += "\ni: " + std::to_string(vi); 
    	for (const vector<int>& vis : G[vi].strip_vis)
            for (int vertex : vis){
                if (q.vertexDegree[vertex] != 4 && !q.boundaryVertices[vertex]){
	                G[vi].label += "\n**" + std::to_string(vertex) + "**"; //" (idx " + to_string(np3dp->vertexDegreeQuad[vertex]) + ")";
                    if (!q.strip_networks.blockedStrips[G[vi].si])
                		G[vi].shape = "hexagon";
                } 
            }
    }

    // write current graph to .dot file
    std::filebuf fb;
    fb.open(dot_filename, std::ios::out);
    std::ostream os(&fb);
    // boost::write_graphviz(os, G, boost::make_label_writer(boost::get(&GraphVertexProps::label, G))); 
    
    // --- write graphviz dot file
	boost::dynamic_properties dp;
	dp.property("color", get(&GraphVertexProps::color, G));
	dp.property("node_id", get(boost::vertex_index, G));
	//dp.property("node_id", get(&GraphVertexProps::si, G));
    dp.property("label", boost::get(&GraphVertexProps::label, G));
    dp.property("shape", get(&GraphVertexProps::shape, G));
    dp.property("label", get(&GraphEdgeProps::replacing_strip, G)); // https://graphviz.org/docs/edges/
    write_graphviz_dp(os, G, dp);
    fb.close();

    // --- convert dot file to png (need to install graphviz for this to work; https://graphviz.org/download/)
	// command : "dot -Kneato -Tpng -Goverlap=scale -Npenwidth=2.0 Dgraph.dot -o _Dgraph.png";
	// More info on command line options: https://graphviz.org/doc/info/command.html (ex. check -klayout)
    string executable = "\"C:\\Program Files\\Graphviz\\bin\\dot.exe\""; // escape double quotes are used to overcome the white space in Program Files : https://stackoverflow.com/questions/17597752/executing-batch-file-from-c-with-spaces-in-the-path
	string command = " -Kneato -Tpng -Goverlap=scale -Npenwidth=3.0 " + dot_filename + " -o " + png_filename; // " -Kneato -Tpng " + 

	system((executable + command).c_str());
}


void Strips::find_blocked_strips(Np3dpContext* np3dp, const vector<AlignedVertices>& aligned_singularities)
{
    QuadMesh& q = np3dp->quad_mesh;

    q.strip_networks.blockedStrips.setZero(q.strip_networks.StoD.rows()); // set all strips as not-blocked

    for (const AlignedVertices& a : aligned_singularities)
	    for (int i=0; i<a.aligned.size(); ++i)
		    if (a.aligned[i])
			    for (int si : a.all_strip_paths[i]) {
                    if (si >= 0){
                        q.strip_networks.blockedStrips[si] = true;
                        GraphVertex gi = q.strip_networks.StoG[si];
                        q.strip_networks.cDgraph[gi].color = "red";
                        q.strip_networks.cDgraph[gi].shape = "rect";
                    }
			    }
    cout << "Marked " << q.strip_networks.blockedStrips.sum() << " strips as blocked." << endl;
}


double Strips::face_squareness(int fi, const MatrixXd& V, const MatrixXi& F, const VectorXi& D) // does not consider direction
{
    double d0_sum = 0; double d1_sum = 0;
    for (int k=0; k<D[fi]; ++k)
    {
	    int v1 = F(fi, k);
        int v2 = F(fi, (k+1)%D[fi]);
        double d = (V.row(v1) - V.row(v2)).norm();
        if (k%2 == 0) d0_sum += d;
        else          d1_sum += d;
    }
    double sq = d0_sum / d1_sum;
    if (sq > 1) sq = 1.0 / sq; // make sq value be between [0.0 and 1.0]
    sq = 1.0 - sq; // 0 if perfectly square, 1 if not square at all
    return sq;
}


double Strips::face_squareness_with_direction(int fi, const MatrixXd& V, const MatrixXi& F, const MatrixXi& FE, const VectorXi& D, const VectorXi& Emap, int dir) // does not consider direction
{
    double d0_sum = 0; double d1_sum = 0;
    for (int k=0; k<D[fi]; ++k)
    {
        int ei = FE(fi, k);
    	int v1 = F(fi, k);
        int v2 = F(fi, (k+1)%D[fi]);
        double d = (V.row(v1) - V.row(v2)).norm();

        if (Emap[ei] == dir)
	        d1_sum += d;
        else
            d0_sum += d;
    }
    double sq = d0_sum / d1_sum;
    sq = sq - 1; // 0 if perfectly square, sq<0: dir too large, sq>0: other dir too large
    return sq;
}


double Strips::face_angles(int fi, const MatrixXd& V, const MatrixXi& F, const VectorXi& D)
{
    double a_sum = 0;
    for (int k = 0; k < D[fi]; ++k)
    {
        int v0 = F(fi, (k-1 + D[fi])%D[fi]);
        int v1 = F(fi, k);
        int v2 = F(fi, (k + 1) % D[fi]);

        RowVector3d vec1 = (V.row(v0) - V.row(v1)).normalized();
        RowVector3d vec2 = (V.row(v2) - V.row(v1)).normalized();

        a_sum += abs(vec1.dot(vec2));
    }
    return a_sum;
}


int Strips::decide_on_next_step(Np3dpContext* np3dp)
{
    QuadMesh& q = np3dp->quad_mesh;
    StripNetworks& sn = q.strip_networks;

    if (q.selection.selected_vis.size() < 2) return -1;

    // --- find intermediate strips between two selected vertices (np3dp->selected_vi1, np3dp->selected_vi2)
    int vi = q.selection.selected_vis[0];
    int vj = q.selection.selected_vis[1];
    vector<bool> aligned;
    vector<vector<int>> alignment_mesh_paths;
    vector<vector<int>> vis_sequences;
    Strips::find_strips_sequences_between_two_vertices(np3dp, vi, vj, sn.counterstrip_sequencies, vis_sequences, aligned, sn.strip_sequencies, alignment_mesh_paths);

    if (std::all_of(aligned.begin(), aligned.end(), [](bool item)->bool {return item == true; })){cout << "Vertices already aligned in all blue directions" << endl; return -1; }// if already aligned, do nothing


    // --- get max thresholds
    int nF = q.F.rows();


    // --- concatenate all strips in one vector
    vector<int> intermediate_strip_indices;
    for (int si : sn.counterstrip_sequencies[sn.selected_sequence_index])
        intermediate_strip_indices.push_back(si);
    for (int si : sn.strip_sequencies[sn.selected_sequence_index])
        intermediate_strip_indices.push_back(si);
    cout << "Chosing strip to use for next step out of " << intermediate_strip_indices.size() << " strips." << endl;

    vector<double> scores(intermediate_strip_indices.size());

    // --- score based on number of quads
    for (int i = 0; i < intermediate_strip_indices.size(); ++i) {
        int si = intermediate_strip_indices[i];

    	if(!sn.blockedStrips[si] || sn.StoD[si] == 0) {
            const VectorXi& strip = sn.strip_index_to_properties(si).strip;
            int strip_size = strip.sum();
    		scores[i] = static_cast<double>(strip_size) / static_cast<double>(nF) ;
    	}
        else
        {
            scores[i] = 1000.0; // avoid that this is ever picked
        }
    }

	// --- score based on how 'nice' the strip is
    for (int i = 0; i < intermediate_strip_indices.size(); ++i) {
        int si = intermediate_strip_indices[i];
        const VectorXi& strip = sn.strip_index_to_properties(si).strip;
        int sum = strip.sum();

        if (sn.StoD[si] == 0) // --- strips
        {
            // evaluate strips based on squareness (because collapsing mostly alters squareness)
            for (int fi = 0; fi < q.F.rows(); ++fi) {
                if (strip[fi] && q.D[fi] == 4) { // only consider quads for this 
                    double sq = face_squareness(fi, q.V, q.F, q.D);
                    scores[i] += 0.1 * sq / static_cast<double>(sum);
                }
            }
        }
        else // --- counter strips
        {
            // evaluate counter-strips based on squares angles (because rewiring mostly alters angles)  
            for (int fi = 0; fi < q.F.rows(); ++fi) {
                if (strip[fi] && q.D[fi] == 4) { // only consider quads for this 
                    double a = face_angles(fi, q.V, q.F, q.D);
                    scores[i] += 0.1 * a / static_cast<double>(sum);
                }
            }
        }
    }

    // --- select strip with smallest score
    const int minElementIndex = std::min_element(scores.begin(), scores.end()) - scores.begin();

    int selected_strip = -1;
    if (scores[minElementIndex] < 3.0)
        selected_strip = intermediate_strip_indices[minElementIndex];

    if (selected_strip >= 0) cout << "Found strip si = " << selected_strip << ", in dir = " << sn.StoD[selected_strip] << " for next step." << endl;

	return selected_strip;
}


//////////////////////////////
// Graph
//////////////////////////////

namespace GraphHelpers
{
    vector<vector<GraphVertex>> get_graph_adjacency_lists(const StripGraph& G){
        const int nV = boost::num_vertices(G);
        vector<vector<GraphVertex>> adj(nV);
        for (int vertex = 0; vertex < nV; ++vertex) {
            for (int v : boost::make_iterator_range(boost::adjacent_vertices(vertex, G)))
                adj[vertex].push_back(v);
        }
        return adj;
    }

    void get_all_paths_util(int nV, GraphVertex u, GraphVertex d, vector<bool>& visited, vector<vector<GraphVertex>>& all_paths, int& current_node_index, int& current_path_index, const vector<vector<GraphVertex>>& adj)
    {
        // Mark the current node and store it in path[]
        visited[u] = true;

        if (all_paths.size() <= current_path_index) { // create path if it doesn't already exist
            all_paths.push_back(vector<GraphVertex>());
            all_paths.back().assign(nV, -1);

            // copy part of previous path up to current_node_index
            int nP = all_paths.size();
            if (nP > 1) // if there is a previously completed path
	            for (int i=0; i< current_node_index; ++i)
	                all_paths.back()[i] = all_paths[nP-2][i];
        }
        vector<GraphVertex>& path = all_paths[current_path_index];
    	path[current_node_index] = u;
        current_node_index++;

        // If current vertex is same as destination then current path has been has been completed, switch to next path
        if (u == d) {
            path.resize(current_node_index); // cut path down to its actual length
            ++current_path_index;
            /*// printout completed path
            cout << "Completed path : ";
            for (int i = 0; i < path.size(); i++)
                cout << path[i] << " ";
            cout << endl;*/
        }
        else // If current vertex is not destination then continue
        {
            // Recur for all the vertices adjacent to current vertex
            for (int v : adj[u])
                if (!visited[v])
                    get_all_paths_util(nV, v, d, visited, all_paths, current_node_index, current_path_index, adj);
        }

        // Retract backwards: Remove current vertex from path[] and mark it as unvisited
        current_node_index--;
        visited[u] = false;
    }
}

int Graph::all_paths_between_two_nodes(const StripGraph& G, GraphVertex source, GraphVertex destination, vector<vector<GraphVertex>>& all_paths)  // Reference : https://www.geeksforgeeks.org/find-paths-given-source-destination/
{
    if (source == destination) {
        all_paths.push_back(vector<GraphVertex>());
    	//all_paths.back().push_back(source);
    }
    else
    {
        // --- create adjacency list from G
        const int nV = boost::num_vertices(G);
        vector<vector<GraphVertex>> adj = GraphHelpers::get_graph_adjacency_lists(G);

        // Mark all the vertices as not visited
        vector<bool> visited; visited.assign(nV, false);

        int current_node_index = 0;
        int current_path_index = 0;

        // Call the recursive helper function to print all paths
        GraphHelpers::get_all_paths_util(nV, source, destination, visited, all_paths, current_node_index, current_path_index, adj);

        // discard last path if it hasn't reached the destination
        if (all_paths.back().back() != destination)
            all_paths.pop_back();
    }

    return all_paths.size();
}


bool Graph::shortest_path(const StripGraph& G, GraphVertex src, const vector<GraphVertex>& targets, vector<int>& dist, vector<int>& pred)
{ // implementing this here instead of using the default boost function, because that function creates a conflict with some graphviz boost function, and the solution won't compile

    // create adjacency list from G
    const int nV = boost::num_vertices(G);
    vector<vector<GraphVertex>> adj = GraphHelpers::get_graph_adjacency_lists(G);

    auto vertex_found = [](const vector<GraphVertex>& targets, GraphVertex vi)->bool {
        for (const auto& v : targets)
            if (v == vi)
                return true;
        return false;
    };

    // Fill in distance and predecessor lists
    //more details here: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/
    typedef pair<int, int> iPair;
    priority_queue<iPair, vector<iPair>, greater<iPair>> pq;

    dist.resize(nV, INF);     // Create a vector for distances and initialize all distances as infinite (INF)
    pred.resize(nV, -1);

    pq.push(make_pair(0, src)); // Insert source itself in priority queue and initialize its distance as 0.
    dist[src] = 0;

    if (vertex_found(targets, src)) // if the src is already in the targets, shortest path has been found
        return true;

    /* Looping till priority queue becomes empty (or targets are found) */
    while (!pq.empty()) {
        // The first vertex in pair is the minimum distance vertex, extract it from priority queue. vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices sorted distance (distance must be first item in pair)
        int u = pq.top().second;
        pq.pop();

        for (int v : adj[u]) { // iterate through all neighbors of u
            if (dist[v] > dist[u] + 1) {
                pred[v] = u;
                dist[v] = dist[u] + 1;
                pq.push(make_pair(dist[v], v));

                if (vertex_found(targets, v))
                    return true;
            }
        }
    }
    return false;
}

