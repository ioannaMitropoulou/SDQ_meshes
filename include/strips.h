//
// Created by Ioanna Mitropoulou on 11.05.22.
//

#pragma once
#ifndef NP_3DP_STRIP_H
#define NP_3DP_STRIP_H

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <Eigen/Core>
#include "Np3dpContext.hpp"


using namespace Eigen;
using namespace std;

struct ShortestPathsResults {
	ShortestPathsResults(GraphVertex node_start, GraphVertex node_end, int si_start, int si_end, int dist, const deque<GraphVertex>& path) : node_start(node_start), node_end(node_end), strip_start(si_start), strip_end(si_end), dist(dist), path(path) {}
	GraphVertex node_start, node_end;
	int strip_start, strip_end;
	int dist;
	deque<GraphVertex> path;
};


struct GraphVisitor : boost::default_dijkstra_visitor {
	using base = boost::default_dijkstra_visitor;
	struct done {};
	GraphVisitor(GraphVertex vd, size_t& visited) : destination(vd), visited(visited) {}

	void finish_vertex(GraphVertex v, StripGraph const& G) {
		++visited;
		if (v == destination)
			throw done{};
		base::finish_vertex(v, G);
	}
private:
	GraphVertex destination;
	size_t& visited;
};


namespace Strips
{
	void get_all_strips(const MatrixXd& V, const MatrixXi& F, const VectorXi& D, const VectorXi& boundaryFaces, const VectorXi& Emap,
		MatrixXi& strips, VectorXi& strips_directions, vector<array<vector<int>, 2>>& strips_vis, vector<bool>& strips_are_closed, vector<vector<int>>& eis,
		vector<vector<int>>& fis, MatrixXi& FtoS, vector<vector<set<int>>>& VtoS, const vector<int>& directions = {0,1}); // by default find strips in both directions

	void create_strip_graphs(Np3dpContext* np3dp);

	void find_strips_sequences_between_two_vertices(Np3dpContext* np3dp, int vi, int vj, vector<vector<int>>& counterstrip_sequences, vector<vector<int>>& vis_sequences, vector<bool>& aligned, vector<vector<int>>& strip_sequences, vector<vector<int>>& alignment_mesh_paths);

	void find_tracked_alignments(Np3dpContext* np3dp, vector<AlignedVertices>& aligned_vis);

	int get_all_counterstrip_sequences_between_two_vertices(const Np3dpContext* np3dp, int vi, int vj, int direction, const VectorXi& Emap, vector<vector<int>>& counterstrip_sequences, vector<vector<int>>& vis_sequences);
	// direction: 0 to get green counter-strips
	// C_intersection_vis: for each path of counterstrips save the vertex where the blue direction from vi intersected with the green direction from vj (i.e. the 3rd corner of the rectangle between vi and vj)

	vector<bool> two_vertices_are_aligned(const Np3dpContext* np3dp, int vi, int vj, int direction, const VectorXi& Emap, const vector<vector<int>>& counterstrip_sequences, const vector<vector<int>>& vis_sequences, vector<vector<int>>& mesh_paths);
	// returns N booleans, one for each of the possible alignment paths between the vertices
	// all_counterstrip_sequences : N paths of counter-strip graph indices that connect the two vertices
	// all_mesh_paths: N paths of quad mesh vertices that connect the two vertices (if aligned[i] == false, then all_mesh_paths[i].empty() == true)

	void write_graphviz(const Np3dpContext* np3dp, const std::string& dot_filename, const std::string& png_filename, StripGraph& G);

	deque<GraphVertex> shortest_strips_path_between_two_mesh_vertices(const Np3dpContext* np3dp, int v1, int v2, int direction, int side_id);
	// path holds *Graph-vertices* of either the Dgraph (if direction==0) or the cDgraph (if direction==1) 

	bool get_strip_from_fi(const VectorXi& D, const MatrixXi& F , const VectorXi& boundaryFaces, const VectorXi& Emap, int fi, 
		                   int strip_direction, VectorXi& strip, array<vector<int>, 2>& strip_vis, bool& strip_is_closed, vector<int>& strip_eis, 
						   vector<int>& strip_fis, bool update_f0_if_boundary = true);

	void find_blocked_strips(Np3dpContext* np3dp, const vector<AlignedVertices>& aligned_singularities);

	double face_squareness(int fi, const MatrixXd& V, const MatrixXi& F, const VectorXi& D);
	double face_squareness_with_direction(int fi, const MatrixXd& V, const MatrixXi& F, const MatrixXi& FE, const VectorXi& D, const VectorXi& Emap, int dir);
	double face_angles(int fi, const MatrixXd& V, const MatrixXi& F, const VectorXi& D);

	int decide_on_next_step(Np3dpContext* np3dp); // we assume that we always align along the dominant direction (The user can switch to other side to align in the other direction)
};


namespace Graph{
	bool shortest_path(const StripGraph& G, GraphVertex src, const vector<GraphVertex>& targets, vector<int>& dist, vector<int>& pred);
	int all_paths_between_two_nodes(const StripGraph& G, GraphVertex source, GraphVertex destination, vector<vector<GraphVertex>>& all_paths); // graph, source, destination
}


#endif
