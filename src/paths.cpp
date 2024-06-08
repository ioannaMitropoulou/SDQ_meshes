#include "paths.h"

#include "helpers.h"


//////////////////////////////////////
/// --- Segment
//////////////////////////////////////

void Paths::Segment::switch_v0_v1() {
    auto cp_v0 = v0;
    auto cp_parentEdges_0 = parentEdges[0];

    v0 = v1;
    parentEdges[0] = parentEdges[1];

    v1 = cp_v0;
    parentEdges[1] = cp_parentEdges_0;
}


//////////////////////////////////////
/// --- Path
//////////////////////////////////////

void  Paths::Path::reverse() {
    std::reverse(segments.begin(), segments.end());
    for (auto& seg : segments)
        seg->switch_v0_v1();
}


vector<int> Paths::Path::get_intersecting_faces() const {
    vector<int> fs;
    for (auto& seg : segments)
        //if (!seg->is_on_sheet_crossing) // do not allow segments that are on a sheet crossing (their face assignment is only by convention)
            fs.push_back(seg->parentFace);
    return fs;
}


vector<int> Paths::Path::get_intersecting_edges() const {
    vector<int> eis;
    for (auto& seg : segments)
        for (int m = 0; m < 2; ++m)
            if (!Helpers::in_vector(eis, seg->parentEdges[m]))
                eis.push_back(seg->parentEdges[m]);
    return eis;
}


shared_ptr<Paths::Segment> Paths::Path::get_segment_on_face(int fi) const {
    for (auto seg : segments)
        if (seg->parentFace == fi)
            //if (!seg->is_on_sheet_crossing) // do not allow segments that are on a sheet crossing (their face assignment is only by convention)
                return seg;

    cerr << "Face fi=" << fi << "is not intersected by path with faces : ";
    Helpers::print_vector(get_intersecting_faces());
    throw;
}


double Paths::Path::avg_dist_from_ground() const {
    double d = 0;
    for (auto seg : segments)
        d += seg->center()[2];
    return d / double(segments.size());
}


double Paths::Path::total_length() const {
    double len = 0;
    for (auto seg : segments)
        len += (seg->v0 - seg->v1).norm();
    return len;
}


void Paths::Path::serialize(int side_id, int piece_id, int path_index, const std::string& serialize_folder) const {
    const std::string filename_start = DATA_PATH + serialize_folder + "Side_" + std::to_string(side_id) + "_Piece_" + std::to_string(piece_id) + "Paths.Path_" + std::to_string(path_index);

    int nSegments = segments.size();
    igl::serialize(nSegments, filename_start + ".segments_size.igl");

    const std::string fileName = filename_start + ".segments.igl";
    std::ofstream f(fileName);
    for (int i = 0; i < nSegments; ++i) {
        shared_ptr<Segment> s = segments[i];
        f << i << " " << s->v0[0] << " " << s->v0[1] << " " << s->v0[2] << " " << s->v1[0] << " " << s->v1[1] << " " << s->v1[2] << " "
            << s->parentEdges[0]  << " " << s->parentEdges[1] << " " << s->parentFace << " "
            << s->forward_dir[0]  << " " << s->forward_dir[1] << " " << s->forward_dir[2] << " "
            << s->normal[0]       << " " << s->normal[1]      << " " << s->normal[2] << 
            endl;
    }
    f.close();
    // assert(!f.fail());
}


void Paths::Path::deserialize(int side_id, int piece_id, int path_index, const std::string& serialize_folder) {
    const std::string filename_start = DATA_PATH + serialize_folder + "Side_" + std::to_string(side_id) + "_Piece_" + std::to_string(piece_id) + "Paths.Path_" + std::to_string(path_index);

    int nSegments;
    igl::deserialize(nSegments, filename_start + ".segments_size.igl");

    //segments.resize(nSegments); // initialize empty
    for (int i = 0; i < nSegments; ++i) {
        Segment s;
        segments.push_back(make_shared<Segment>(s));
    }

    const std::string fileName = filename_start + ".segments.igl";
    std::string line;

    ifstream myfile(fileName);
    if (myfile.is_open()) {
        int i=0;
        while (getline(myfile, line))
        {
            std::istringstream in(line);
            int parentEdges0, parentEdges1, parentFace;
            float v0x, v0y, v0z, v1x, v1y, v1z, g0, g1, g2, n0, n1, n2;
            in >> i >> v0x >> v0y >> v0z >> v1x >> v1y >> v1z >> parentEdges0 >> parentEdges1 >> parentFace >> g0 >> g1 >> g2 >> n0 >> n1 >> n2;
            segments[i]->v0 = Vector3d(v0x, v0y, v0z);
            segments[i]->v1 = Vector3d(v1x, v1y, v1z);
            segments[i]->parentEdges = Vector2i(parentEdges0, parentEdges1);
            segments[i]->parentFace = parentFace;
            segments[i]->forward_dir = Vector3d(g0, g1, g2);
            segments[i]->normal = Vector3d(n0, n1, n2);
        }
        myfile.close();
        if (i != 0 && i + 1 != nSegments) { cerr << "Attention! Wrong number of segments in path found while deserializing data" << endl; throw; } 
    }
    else {
        cerr << "Could not deserialize file : " << filename_start + ".segments.igl" << endl; throw;
    }

    if (segments.size()>=2) // check if path is closed
    {
    	if (segments[0]->parentEdges[0] == segments.back()->parentEdges[1])
            is_closed = true;
    }
}


int Paths::Path::closest_segment_index_to_pt(Vector3d pt, double& min_dist) const
{
    int nSeg = segments.size();
    vector<double> ds(nSeg);
    for (int seg_i = 0; seg_i < nSeg; ++seg_i)
    {
	    shared_ptr<Paths::Segment> c_seg = segments[seg_i];
        Vector3d cpt = Geometry::project_pt_on_line_segment(c_seg->v0, c_seg->v1, pt);
        ds[seg_i] = (cpt - pt).norm();
    }
	const int closest_seg_i = std::distance(ds.begin(), std::min_element(ds.begin(), ds.end())); // closest segment index

    min_dist = *std::min_element(ds.begin(), ds.end()); // min distance
    return closest_seg_i;
}



//////////////////////////////////////
/// --- PathCollection
//////////////////////////////////////

int Paths::PathCollection::add(const Path& p) {
    paths.push_back(p);
    return paths.size() - 1;
}


void Paths::PathCollection::write_to_file(const std::string& fileName) {
    std::ofstream f(fileName);

    for (int pi=0; pi<paths.size(); ++pi) {
        f << "is_closed , " << int(paths[pi].is_closed) << endl;
        for (auto seg : paths[pi].segments)
        {
            f << "s , " << seg->v0[0] << " , " << seg->v0[1] << " , " << seg->v0[2] << " , "
                << seg->v1[0] << " , " << seg->v1[1] << " , " << seg->v1[2]  << " , "
        		<< seg->normal[0] << " , " << seg->normal[1] << " , " << seg->normal[2] << " , " 
                << seg->parentFace  << std::endl;
        }
        f << "end , " << std::endl;
    }
    f.close();
    std::cout << "Saved paths collection : " << fileName << std::endl;
}


void Paths::PathCollection::display(igl::opengl::glfw::Viewer& viewer, int& contours_viewer_index, bool clear_existing=true) {
    vector<shared_ptr<Paths::Segment>> segments;

    for (int pi = 0; pi<paths.size(); ++pi) {
        const Path& p = paths[pi];
        for (auto& seg : p.segments)
            segments.push_back(seg);
    }
    int nS = segments.size();
    VectorXd color_attr;
    color_attr.setZero(nS);

    MatrixXd segment_colors;
    segment_colors.setZero(nS, 3);

    MatrixXd P1(nS, 3), P2(nS, 3); // all start (P1) and end(P2) points of the iso-contours segments : (No_of_line_segments x 3)
    for (int i = 0; i < nS; ++i) {
        P1.row(i) = segments[i]->v0;
        P2.row(i) = segments[i]->v1;
    }

    if (contours_viewer_index < 0) // if viewer index has not been initialized.
    {
        viewer.append_mesh(); // viewer.selected_data_index = 1
        contours_viewer_index = viewer.selected_data_index;
        viewer.data_list[contours_viewer_index].line_width = 1.0;
        viewer.data_list[contours_viewer_index].add_edges(P1, P2, segment_colors); //  colors : RowVector3d(0.0, 0.0, 0.0)

        //viewer.data_list[contours_viewer_index].point_size = 4.0;
    	//viewer.data_list[contours_viewer_index].add_points(P1, RowVector3d(0.0, 0.0, 0.0));
        //viewer.data_list[contours_viewer_index].add_points(P2, RowVector3d(0.0, 0.0, 0.0));

        viewer.selected_data_index = 0; // return to default value
    }
    else // if a contours drawing already exists
    {
        if (clear_existing) viewer.data_list[contours_viewer_index].clear(); // clearing again? (You don't really need to do this here, it has been cleared in np-3dp)
        viewer.data_list[contours_viewer_index].add_edges(P1, P2, segment_colors); //colors    Eigen::RowVector3d(0.0, 0.0, 0.0)

        //viewer.data_list[contours_viewer_index].point_size = 4.0;
        //viewer.data_list[contours_viewer_index].add_points(P1, RowVector3d(0.0, 0.0, 0.0));
        //viewer.data_list[contours_viewer_index].add_points(P2, RowVector3d(0.0, 0.0, 0.0));
    }
}


void Paths::PathCollection::serialize(int side_id, int piece_id, const std::string& serialize_folder) const {
    int nPaths = paths.size();
    igl::serialize(nPaths, DATA_PATH + serialize_folder +  + "Side_" + std::to_string(side_id) + "_Piece_" + std::to_string(piece_id) + "Paths.PathCollection.nPaths.igl");
    for (int i = 0; i < paths.size(); ++i) {
        auto& p = paths[i];
        p.serialize(side_id, piece_id, i, serialize_folder);
    }
}


void Paths::PathCollection::deserialize(int side_id, int piece_id, const std::string& serialize_folder) {
    paths.clear();
    int nPaths;
    igl::deserialize(nPaths, DATA_PATH + serialize_folder +  + "Side_" + std::to_string(side_id) + "_Piece_" + std::to_string(piece_id) + "Paths.PathCollection.nPaths.igl");
    for (int i = 0; i < nPaths; ++i) {
        Path p;
        paths.push_back(p);
        paths.back().deserialize(side_id, piece_id, i, serialize_folder);
    }
}


bool Paths::PathCollection::check_edge_links()
{
	for (auto& p : paths) {
		for (int i=1; i<p.segments.size(); ++i) {
            if (p.segments[i-1]->parentEdges[1] != p.segments[i]->parentEdges[0]) {
                cout << " --- Attention! Not matching edge links at PathCollection! \np.segments[i-1].parentFace = " << p.segments[i - 1]->parentFace
            	<< "\np.segments[i].parentFace = " << p.segments[i]->parentFace << "\np.segments[i-1].parentEdges = " << p.segments[i - 1]->parentEdges
            	<< "\np.segments[i].parentEdges = " << p.segments[i]->parentEdges << endl;
                //throw;
                return false;
            }
		}
	}
    return true;
}