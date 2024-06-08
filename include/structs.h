//
//  Created by Ioanna Mitropoulou
//

#pragma once
#ifndef NP_3DP_STRUCTS_H
#define NP_3DP_STRUCTS_H

#include <Eigen/Core>

#include "paths.h"
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;
using namespace Eigen;

// --- Boost graph definitions
struct GraphVertexProps {
    VectorXi strip; // strip vector (#Fx1) with 1/0 for faces that belong /don't belong to the strip
    bool strip_is_closed = false;
    std::array<vector<int>, 2> strip_vis; // one vector of vis for each Strip Network of the strip
    int si = -1; // strip index

    // display properties
    std::string label;
    std::string color = "blue";
    std::string shape = "oval"; // rect, hexagon // see more shapes here: https://graphviz.org/doc/info/shapes.html
	// other attributes: https://graphviz.org/docs/nodes/
};

struct GraphEdgeProps {
	int replacing_strip = -1;
};

using StripGraph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, GraphVertexProps, GraphEdgeProps>;
using GraphVertex = boost::graph_traits<StripGraph>::vertex_descriptor;
using Edge = boost::graph_traits<StripGraph>::edge_descriptor;
using IndexMap = boost::property_map<StripGraph, boost::vertex_index_t>::type;
using vertex_iter = boost::graph_traits<StripGraph>::vertex_iterator;
using GraphTraits = boost::graph_traits<StripGraph>;


// --- Singularity alignments
struct AlignedVertices {
    AlignedVertices(int vi, int vj, int dir, const vector<bool>& aligned, const vector<vector<int>>& all_mesh_paths, const vector<vector<int>>& all_strip_paths):
	direction(dir), aligned(aligned), all_mesh_paths(all_mesh_paths), all_strip_paths(all_strip_paths){
	    vis[0] = vi;
    	vis[1] = vj;
        N = aligned.size();
        if (all_mesh_paths.size() != N || all_strip_paths.size() != N){ cerr << "The vectors aligned, all_mesh_paths, all_graph_paths should all have the same size" << endl; throw;}
    }
    array<int, 2> vis;
    int direction;
    vector<bool> aligned;
    vector<vector<int>> all_mesh_paths;
    vector<vector<int>> all_strip_paths; // path of counterstrips
    int N; // number of alignment paths
    void append_visualization_data(const MatrixXd& Vquad, MatrixXd& P1, MatrixXd& P2, vector<int>& vis_out) const; // appends visualization data of current aligned pair on P1, P2, vis  
};


struct SelectionData
{
    int selected_fi = -1; // selected face index
	vector<int> selected_vis; // selected vertex indices
    int n_vis_to_select = 2; // how many vertices to select
    bool vis_4 = false; // if true, then the user can select 4 vis instead of 2
	int selected_ei = -1; // selected edge index
    int selected_si = -1; // selected strip index
    vector<int> selected_boundary_eis;

	int selected_strip_direction = 0;
    VectorXi selected_strip; // (Fx1) 0: face not on strip, 1: face on strip
    bool strip_is_closed = false;
    array<vector<int>, 2> selected_strip_vis; // one vector of vis for each side of the strip
    vector<int> selected_strip_fis; // list of indices of strip faces in order

    void invalidate();
};

struct StripNetworks 
{ // data related to the two strip networks and their graphs
	bool has_strip_networks_data = false;

    StripGraph Dgraph, cDgraph;
    MatrixXi strips; // SxF : Fx1 vector per strips (one each row) with 1s/0s for faces that are/are not in the strip
    MatrixXi FtoS; // Fx2: face to strip (for each strip-direction)
    vector<vector<set<int>>> VtoS; // Vx2: vertex to strip (for each strip-direction):  VtoS[vi][strip_direction][si]
    VectorXi StoG; // Strip to Graph vertex index (Strips x 1) - includes both strips in dir0 and dir1
    VectorXi StoD; // Strip to Direction (Strips x 1) - 0 for dominant strips saved in D, 1 for counter-strips saved in cD
    VectorXi EtoS; // Edge to strip index (has -1 if edge is not on the interior of any strip). An edge is mapped to a strip if it is *inside* that strip, and not on its boundary 
	vector<vector<int>> StoEis;// Strip to Edge indices. For each strip: a vector of ordered edge indices of the edges that are *inside* the strip, i.e. orthogonal to its direction
    vector<vector<int>> StoFis;// Strip to Face indices. For each strip: a vector of ordered face indices

	VectorXi blockedStrips; // (Strips x 1) 1/0 for strips that are/aren't blocked for rewiring - includes both strips in dir0 and dir1

    int selected_sequence_index = 0; // from all (counter)strip sequences between two vertices, this stores which sequence is selected
    vector<vector<int>> strip_sequencies; // all strip sequences between two vertices 
    vector<vector<int>> counterstrip_sequencies; // all couterstrip sequences between two vertices. For each strip sequence we have one counterstrip sequence, strip_sequencies.size() == counterstrip_sequencies.size()
    vector<vector<int>> vis_sequences; // holds the sequence of vertex indices that are visited while following each green sequence from v1 to v2
    VectorXi strip_to_modify; // strip index selected for modification

    const GraphVertexProps& strip_index_to_properties(int si) const;

    void invalidate()
    {
	    has_strip_networks_data = false;        
    	Dgraph = StripGraph();
	    cDgraph = StripGraph();

        strips.setZero(0, 0);
        FtoS.setZero(0, 0);
        VtoS.clear();
        StoG.setZero(0);
        StoD.setZero(0);
        EtoS.setZero(0);
        StoEis.clear();
        StoFis.clear();

        blockedStrips.setZero(0);

        selected_sequence_index = 0;
        strip_sequencies.clear();
        counterstrip_sequencies.clear();
        vis_sequences.clear();
        strip_to_modify.setZero(0);
    }
};


#endif //NP_3DP_STRUCTS_H