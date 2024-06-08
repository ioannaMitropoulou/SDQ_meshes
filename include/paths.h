//
// Created by Ioanna Mitropoulou on 29.11.21.
//

#pragma once
#ifndef NP_3DP_PATHS_H
#define NP_3DP_PATHS_H

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;


namespace Paths {

    struct Segment {
        Segment() = default;

        Vector3d v0, v1; // start and end vertex
        Vector2i parentEdges; // the edges where v0 and v1 lie
        int parentFace = -1; // the face where the segment lies
        Vector3d forward_dir; // the forward direction on this face (orthogonal to the normal and to the path)
        Vector3d normal; // the smooth normal (comes from current and neighboring face)

        void switch_v0_v1();
        inline Vector3d center() const { return 0.5 * (v0 + v1); }
        inline Vector3d vector() const { return (v1 - v0).normalized(); }
    };

    struct Path {
        vector<shared_ptr<Segment>> segments;

        bool is_closed = false;

    	void serialize(int side_id, int piece_id, int path_index, const std::string& serialize_folder) const;
        void deserialize(int side_id, int piece_id, int path_index, const std::string& serialize_folder);

        void reverse();
        vector<int> get_intersecting_faces() const; // for N segments: N fs
        vector<int> get_intersecting_edges() const; // for N segments: N+1 edges for open paths, N edges for closed paths
        shared_ptr<Segment> get_segment_on_face(int fi) const;
        double avg_dist_from_ground() const;
        double total_length() const;
        int closest_segment_index_to_pt(Vector3d pt, double& dist) const;

        int ei_first()     const { return segments[0]->parentEdges[0]; }
        Vector3d v_first() const { return segments[0]->v0; }
        int ei_last() const { return segments.back()->parentEdges[1]; }
        Vector3d v_last() const { return segments.back()->v1; }
    };

    struct PathCollection {
        vector<Path> paths;

    	void serialize(int side_id, int piece_id, const std::string& serialize_folder) const;
        void deserialize(int side_id, int piece_id, const std::string& serialize_folder);

        void insert_paths_at_beginning(const PathCollection& other) { paths.insert(paths.begin(), other.paths.begin(), other.paths.end()); }
        void insert_paths_at_end      (const PathCollection& other) { paths.insert(paths.end(),   other.paths.begin(), other.paths.end()); }

        int add(const Path& p);
        void write_to_file(const std::string& fileName);
        void display(igl::opengl::glfw::Viewer& viewer, int& contours_viewer_index, bool clear_existing);
        bool check_edge_links();

        inline void clear() { paths.clear(); }
        inline void remove(int pi) { paths.erase(paths.begin() + pi); }
    };
}



#endif //NP_3DP_PATHS_H