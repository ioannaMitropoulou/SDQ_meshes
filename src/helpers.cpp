//
// Created by Ioanna Mitropoulou on 29.11.21.
//
#include "helpers.h"

MatrixXd Helpers::read_matrix_Xx3_from_txt(const std::string& fileName)
{
    MatrixXd s;
    string line;
    ifstream myfile(fileName);

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            std::istringstream in(line);
            double x, y, z;
            in >> x >> y >> z;

            int current_row = s.rows();
            s.conservativeResize(current_row + 1, 3);
            s.row(current_row) = Vector3d(x, y, z);
        }
        myfile.close();
    }
    else
    {
        cerr << "Could not read file : " << fileName << endl;
    }
    return s;
}


MatrixXd Helpers::read_matrix_Xx6_from_txt(const std::string& fileName)
{
    MatrixXd s;
    string line;
    ifstream myfile(fileName);

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            std::istringstream in(line);
            double x, y, z, a, b, c;
            in >> x >> y >> z >> a >> b >> c;

            int current_row = s.rows();
            s.conservativeResize(current_row + 1, 6);
            s.row(current_row)  << x, y, z, a, b, c;
        }
        myfile.close();
    }
    else
    {
        cerr << "Could not read file : " << fileName << endl;
    }
    return s;
}


VectorXi Helpers::read_vectorXi_from_txt(const std::string& fileName)
{
    string line;
    ifstream myfile(fileName);

    vector<int> vec;
    VectorXi s;

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {

            std::istringstream in(line);
            int integer;
            in >> integer;

            vec.push_back(integer);
        }
        myfile.close();

        s.setZero(vec.size());
        for (int i = 0; i < vec.size(); ++i)
        {
            s[i] = vec[i];
        }
    }
    else
    {
        cerr << "Could not read file : " << fileName << endl;
    }
    return s;
}


VectorXi Helpers::read_numbered_vectorXi_from_txt(const std::string& fileName)
{
	string line;
	ifstream myfile(fileName);

	VectorXi s;

	if (myfile.is_open())
	{
	    int maxIndex = 0;
	    // First pass to get the maximum index to size our vector correctly
	    while (getline(myfile, line))
	    {
	        std::istringstream in(line);
	        string token;
	        int index, value;

	        getline(in, token, ','); 
	        index = stoi(token); 
	        if(index > maxIndex) {
	            maxIndex = index;
	        }
	    }

	    s.setZero(maxIndex + 1); // assume indices start from 0

	    // Clear EOF and rewind
	    myfile.clear();
	    myfile.seekg(0, ios::beg);

	    // Second pass to read values into the vector
	    while (getline(myfile, line))
	    {
	        std::istringstream in(line);
	        string token;
	        int index, value;

	        getline(in, token, ','); 
	        index = stoi(token); 

	        getline(in, token); 
	        value = stoi(token);

	        s[index] = value; // directly assign value to the index-th position of s
	    }

	    myfile.close();
	}
	else
	{
	    cerr << "Could not read file : " << fileName << endl;
	}
	return s;
}


bool Helpers::read_vf_constraints(const std::string& fileName, VectorXi& fis, MatrixXd& dirs, VectorXd& ws) {
    ifstream data(fileName);
    if (!data.is_open()) {
        cerr << "Could not open input fileName" << fileName << endl;
        return false;
    }

    vector<int> fis_vec;
    vector<Vector3d> dirs_vec;
    vector<double> weights_vec;

    string s;
    while (getline(data, s)) {
        // file format: index, angle, vec_x, vec_y, vec_z
        stringstream line(s);
        int face_index;
        float x, y, z, w;
        line >> face_index >> x >> y >> z >> w;

        fis_vec.push_back(face_index); // append face index to b
        dirs_vec.push_back(Vector3d(x, y, z)); // append vector coordinates to bc
        weights_vec.push_back(w);
    }
    data.close();

    fis.setZero(fis_vec.size());
    dirs.setZero(dirs_vec.size(), 3);
    ws.setZero(weights_vec.size());
    for (int i = 0; i < fis_vec.size(); ++i) {
        int fi = fis_vec[i];
        fis[i] = fi;
        dirs.row(i)= dirs_vec[i];
        ws[i] = weights_vec[i];
    }

    cout << "Loaded directional constraints on " << dirs.rows() << " faces." << endl;

    return true;
}


void Helpers::setup_line_texture(K& texture_R, K& texture_G, K& texture_B, int texture_scale, int lineWidth)
{
    unsigned size = 128;
    texture_B.setConstant(size, size, 255);
    texture_G.setConstant(size, size, 255);
    texture_R.setConstant(size, size, 255);

    int squares_num = texture_scale;
    int squares_size = size / texture_scale;
    int last_square_size = size - ((squares_num - 1) * squares_size);

    if (lineWidth >= min(squares_size, last_square_size))
        lineWidth = int(min(squares_size, last_square_size) * 0.75);

    if (lineWidth == 0)
        cerr << "Attention, lineWidth = 0" << endl;

    for (int square_i = 0; square_i < squares_num; ++square_i) {

        for (int square_j = 0; square_j < squares_num; ++square_j) {

            unsigned start_i = square_i * squares_size;
            unsigned start_j = square_j * squares_size;

            unsigned current_size = square_i == squares_num - 1 ? squares_size : last_square_size;
            unsigned current_size_2 = current_size / 2;

            for (unsigned i = start_i; i < start_i + current_size; ++i)
                for (unsigned j = start_j + current_size_2 - lineWidth; j <= start_j + current_size_2 + lineWidth; ++j) { // orange // green
                    texture_R(i, j) = 0; //255.0;
		            texture_G(i, j) = 0.9 * 255; //215;// 225.0; // 215.0;
		            texture_B(i, j) = 0.2 * 255; //100; // 255.0; // 50.0;
                }
            for (unsigned i = start_j + current_size_2 - lineWidth; i <= start_j + current_size_2 + lineWidth; ++i)
                for (unsigned j = start_i; j < start_i + current_size; ++j) { // black // blue
                    //texture_B(i, j) = texture_G(i, j) = texture_R(i, j) = 0 + 100; //+ 200;
                    texture_R(i, j) = 0; 
                    texture_G(i, j) = 0.2 * 255; 
                    texture_B(i, j) = 0.9 * 255; 
                }
        }
    }
}


double Helpers::remap_unbound(double input_val, double in_from, double in_to, double out_from, double out_to)
{
    double out_range = out_to - out_from;
    double in_range = in_to - in_from;
    double in_val = input_val - in_from;
    double val = (float(in_val) / in_range) * out_range;
    return out_from + val;
}


double Helpers::remap(double input_val, double in_from, double in_to, double out_from, double out_to)
{
    if (input_val < in_from) 
        return out_from;
    else if (input_val > in_to)
         return out_to;
    else
        return remap_unbound(input_val, in_from, in_to, out_from, out_to);
}


void Helpers::write_VV_to_txt(const std::string& fileName, const vector<vector<int>>& VV)
{
    std::ofstream f(fileName);
	for (int vi=0; vi<VV.size(); ++vi){
		for (int vn : VV[vi])
		{
            f << vn << " ";
		}
        f << endl;
	}
    cout << "Saved : " << fileName << endl;
}


void Helpers::write_Emap_to_txt(const string& fileName, const VectorXi& Emap, const MatrixXi& EVquad, bool quiet)
{
    if (EVquad.rows() != Emap.rows()) {cerr << "While writing Emap to txt, wrong size of Emap or EVquad : EVquad.rows() != Emap.rows()" << endl; throw;}
    MatrixXi Emap_export(Emap.rows(), 3);
    for (int ei = 0; ei < Emap.rows(); ++ei)
        Emap_export.row(ei) = RowVector3i(EVquad(ei, 0), EVquad(ei, 1), Emap[ei]);
    Helpers::write_matrix_to_txt(Emap_export, fileName, quiet);
}


void Helpers::remove_directory_contents(const std::string& directory_path)
{

#if defined(_WIN32)
    WIN32_FIND_DATA findFileData;
    HANDLE hFind = FindFirstFile((directory_path + "\\*").c_str(), &findFileData);

    if (hFind == INVALID_HANDLE_VALUE) {
        std::cerr << "Failed to find first file in the directory" << std::endl;
        return;
    }

    do {
        const std::string fileOrDirName = findFileData.cFileName;
        if (fileOrDirName == "." || fileOrDirName == "..") continue;

        std::string fullPath = directory_path + "\\" + fileOrDirName;

        if (findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            // Recursively delete directory contents
            remove_directory_contents(fullPath);
            RemoveDirectory(fullPath.c_str()); // Remove the now-empty directory
        } else {
            DeleteFile(fullPath.c_str()); // Delete the file
        }
    } while (FindNextFile(hFind, &findFileData) != 0);

    FindClose(hFind);

#else
    DIR* dir = opendir(directory_path.c_str());
    struct dirent* entry;

    if (dir == nullptr) {
        std::cerr << "Failed to open directory" << std::endl;
        return;
    }

    while ((entry = readdir(dir)) != nullptr) {
        const std::string fileOrDirName = entry->d_name;
        if (fileOrDirName == "." || fileOrDirName == "..") continue;

        std::string fullPath = directory_path + "/" + fileOrDirName;

        struct stat path_stat;
        stat(fullPath.c_str(), &path_stat);
        if (S_ISDIR(path_stat.st_mode)) {
            // Recursively delete directory contents
            remove_directory_contents(fullPath);
            rmdir(fullPath.c_str()); // Remove the now-empty directory
        } else {
            unlink(fullPath.c_str()); // Delete the file
        }
    }

    closedir(dir);
#endif
}


void Helpers::create_directory_if_not_exists(const std::string& path)
{
    #if defined(_WIN32) || defined(_WIN64)
        // Windows-specific code
        DWORD dwAttrib = GetFileAttributesA(path.c_str());

        if (dwAttrib == INVALID_FILE_ATTRIBUTES || !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)) {
            // Directory does not exist.
            if (!CreateDirectoryA(path.c_str(), NULL)) {
                std::cerr << "Error creating directory: " << path << std::endl;
            } else {
                std::cout << "Created directory: " << path << std::endl;
            }
        } else {
            // Directory exists.
            std::cout << "Directory already exists: " << path << std::endl;
        }
    #else
        // POSIX (Linux, macOS, etc.) specific code
        struct stat info;

        if (stat(path.c_str(), &info) != 0) {
            // Directory does not exist.
            if (mkdir(path.c_str(), 0755) == -1) {
                std::cerr << "Error creating directory: " << path << std::endl;
            } else {
                std::cout << "Created directory: " << path << std::endl;
            }
        } else if (info.st_mode & S_IFDIR) {
            // Directory exists.
            std::cout << "Directory already exists: " << path << std::endl;
        } else {
            // Path exists but is not a directory.
            std::cerr << "Path exists but is not a directory: " << path << std::endl;
        }
    #endif
}




///////////////////////////////////
/// Geometry
///////////////////////////////////

double Geometry::get_min_distance_from_pts(const Vector3d& pt, const vector<Vector3d>& cloud)
{
    double min_d = 1e10;
    for (const Vector3d& mark : cloud) {
        min_d = std::min((mark - pt).squaredNorm(), min_d);
    }
    return pow(min_d, 0.5);
}


Vector3d Geometry::project_pt_on_line(const Vector3d& A, const Vector3d& B, const Vector3d& P) {
    //returns the closest point to P on the line segment from A to B
    Vector3d AB = B - A;

    // Consider the line extending the segment, parameterized as A + t (B - A)
    // We find projection of point p onto the line. It falls where t = [(p-A) . (B-A)] / |B-A|^2
    Vector3d AP = P - A;
    double t = AP.dot(AB) / AB.squaredNorm();
    return A + t * AB;
}


bool Geometry::projected_point_inside_segment(const Vector3d& p, const Vector3d& q, const Vector3d& r)
{ // returns true if q is on the segment p-r.
    if (q[0] <= max(p[0], r[0]) && q[0] >= min(p[0], r[0])
        && q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1])
        && q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2]))
        return true;
    return false;
}


Vector3d Geometry::project_pt_on_plane(const Vector3d& N, const Vector3d& origin, const Vector3d& pt)
{
    return pt - (pt - origin).dot(N) * N;
}


pair<Vector3d, bool> Geometry::intersection_line_plane(const Vector3d& a, const Vector3d& b, const Vector3d& o, const Vector3d& n, double tol) {
    // a,b: line endpoints. o,n: plane origin and normal (unit vector)
    // reference : https://github.com/compas-dev/compas/blob/d91265e7147ddeb1abd25a10e416acee820dcf73/src/compas/geometry/intersections/intersections.py#L232
    Vector3d ab = b - a;
    const double cosa = n.dot(ab);

    if (abs(cosa) <= tol)
        /*if the dot product (cosine of the angle between segment and plane)
        is close to zero the line and the normal are almost perpendicular, hence there is no intersection */
        return make_pair(Vector3d(0, 0, 0), false);

    /*based on the ratio = -dot_vectors(n, ab) / dot_vectors(n, oa)
    there are three scenarios
    1) 0.0 < ratio < 1.0: the intersection is between a and b
    2) ratio < 0.0: the intersection is on the other side of a
    3) ratio > 1.0: the intersection is on the other side of b*/

    const Vector3d oa = o - a;
    const double ratio = n.dot(oa) / cosa;
    ab *= ratio;
    return make_pair(a + ab, true);
}


pair<Vector3d, bool> Geometry::intersection_line_line(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d, double tol = 1e-3) {

    auto cross = [](const Vector3d& a, const Vector3d& b)->Vector3d {
        // When I run the project in Release mode, the Eigen cross function throws a linker error, I can't figure out why.
        // So for now I'm just implementing it here
        // TODO: investigate further.
        Vector3d cross;
        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[2] * b[0] - a[0] * b[2];
        cross[2] = a[0] * b[1] - a[1] * b[0];
        return cross;
    };

	// returns true if the line segments have intersections _inside_ the segments
    const Vector3d ab = b - a;
    const Vector3d cd = d - c;
    //    const Vector3d n = ab.cross(cd);
    //    const Vector3d n1 = ab.cross(n).normalized();
    //    const Vector3d n2 = cd.cross(n).normalized();

    Vector3d n = cross(ab, cd);
    Vector3d n1 = cross(ab, n).normalized(); // plane parallel to the ab vector
    Vector3d n2 = cross(cd, n).normalized(); // plane parallel to the cd vector

    pair<Vector3d, bool> r1 = intersection_line_plane(a, b, c, n2, tol);
    pair<Vector3d, bool> r2 = intersection_line_plane(c, d, a, n1, tol);

    if (r1.second && r2.second) // if both intersect
        if ((r1.first - r2.first).squaredNorm() < tol) // if the two pts are identical
            if (projected_point_inside_segment(a, r1.first, b) && projected_point_inside_segment(c, r1.first, d))
                return make_pair(r1.first, true);
    return make_pair(Vector3d(0, 0, 0), false);
}


pair<Vector3d, bool> Geometry::projected_intersection_line_line(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d,
    const Vector3d& plane_N, const Vector3d& plane_origin)
{
    // first project two lines on the same plane (the average of the normals of the two planes where they currently belong)
    Vector3d a_proj = project_pt_on_plane(plane_N, plane_origin, a);
    Vector3d b_proj = project_pt_on_plane(plane_N, plane_origin, b);
    Vector3d c_proj = project_pt_on_plane(plane_N, plane_origin, c);
    Vector3d d_proj = project_pt_on_plane(plane_N, plane_origin, d);

    // then check for intersection
    return intersection_line_line(a_proj, b_proj, c_proj, d_proj);
}


Vector3d Geometry::project_pt_on_line_segment(const Vector3d& A, const Vector3d& B, const Vector3d& P) {
    //returns q the closest point to p on the line segment from A to B
    // https://www.alecjacobson.com/weblog/?p=1486
    Vector3d AB = B - A;
    double AB_sq = AB.squaredNorm();

    if (AB_sq == 0) { // A and B are the same point
        return A;
    }
    else {
        // Consider the line extending the segment, parameterized as A + t (B - A)
        // We find projection of point p onto the line. It falls where t = [(p-A) . (B-A)] / |B-A|^2
        Vector3d AP = P - A;
        double t = AP.dot(AB) / AB_sq;

        if (t < 0) // "Before" A on the line, just return A
            return A;
        else if (t > 1.0) // "After" B on the line, just return B
            return B;
        else //projection lines "inbetween" A and B on the line
            return A + t * AB;
    }
}


Matrix4d Geometry::rotation_matrix_from_vectors(const Vector3d& vec1, const Vector3d& vec2)
{
    auto cross = [](const Vector3d& a, const Vector3d& b)->Vector3d {
        // When I run the project in Release mode, the Eigen cross function throws a linker error, I can't figure out why.
        // So for now I'm just implementing it here
        // TODO: investigate further.
        Vector3d cross;
        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[2] * b[0] - a[0] * b[2];
        cross[2] = a[0] * b[1] - a[1] * b[0];
        return cross;
    };

    Vector3d v1 = vec1.normalized();
    Vector3d v2 = vec2.normalized();

    if ((v1 - v2).squaredNorm() < 1e-4) {// if the vectors are identical, then no rotation is needed, just return identity
        Matrix4d Id; Id.setIdentity();
        return Id;
    }

    // get orthonormal frame from first vector
    Vector3d y1, x1;
	{
        Vector3d arbitrary_other_vec(1.0, 0.0, 0.0);
        if (abs((arbitrary_other_vec - v1).squaredNorm() - 1) < 1e-6 || abs((arbitrary_other_vec + v1).squaredNorm() - 1) < 1e-6) { arbitrary_other_vec = Vector3d(0.0, 0.0, 1.0).normalized(); }
    	y1 = cross(v1, arbitrary_other_vec).normalized();
    	x1 = cross(v1, y1).normalized();
	}
    //cout << "v1 : " << v1 << endl;
    //cout << "y1 : " << y1 << endl;
    //cout << "x1 : " << x1 << endl;

    // get orthonormal frame from second vector
    Vector3d y2, x2;
	{
		Vector3d arbitrary_other_vec(1.0, 0.0, 0.0);
    	if (abs((arbitrary_other_vec - v2).squaredNorm() - 1) < 1e-6 || abs((arbitrary_other_vec + v2).squaredNorm() - 1) < 1e-6) { arbitrary_other_vec = Vector3d(0.0, 0.0, 1.0).normalized(); }
    	y2 = cross(v2, arbitrary_other_vec).normalized();
    	x2 = cross(v2, y2).normalized();
	}
    //cout << "v2 : " << v2 << endl;
    //cout << "y2 : " << y2 << endl;
    //cout << "x2 : " << x2 << endl;

    // get transformations
    Matrix4d R1; R1.setIdentity(); // transforms from world to (v1,y1,z1)
    R1.block(0,0,3,1) = x1;
    R1.block(0, 1, 3, 1) = y1;
    R1.block(0, 2, 3, 1) = v1;

    Matrix4d R2; R2.setIdentity();  // transforms from world to (v2,y2,z2)
    R2.block(0, 0, 3, 1) = x2;
    R2.block(0, 1, 3, 1) = y2;
    R2.block(0, 2, 3, 1) = v2;

	/*if (!(R1.transpose() * R1).isIdentity())
	{
		cerr << "R1.transpose() * R1 is not identity!" << endl;
        cout << "R1 : " << endl << R1 << endl;
        cout << "R1.transpose() * R1 : " << endl << R1.transpose() * R1 << endl;
        throw;
	}
    if (!(R2.transpose() * R2).isIdentity())
    {
	    cerr << "R2.transpose() * R2 is not identity!" << endl;
        cout << "R2 : " << endl << R2 << endl;
        cout << "R2.transpose() * R2 : " << endl << R2.transpose() * R2 << endl;
        throw;
    }*/

    return R2 * R1.transpose();
}


double Geometry::distance_of_pt_from_plane(const RowVector3d& vertex, const RowVector3d& N, const RowVector3d& C)
{
	return std::abs((vertex - C).dot(N)) / N.norm();
}

bool Geometry::PCA(const MatrixXd& P, VectorXd& eigenvalues, MatrixXd& eigenvectors)
{
    MatrixXd centered = P.rowwise() - P.colwise().mean();
    MatrixXd cov = centered.transpose() * centered;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
    if (eigensolver.info() != Eigen::Success) { cerr << "Eigen decomposition failed" << endl; return false; }
    eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();

    // the eigenvalues do not correspond with the extents of the geometry. (I thought that they should..). So below I just find myself the extents of geometry on each eigenvector direction
    /* {
        if (P.cols() != eigenvectors.cols()) { cerr << "The eigenvectors do not have the same number of columns as P" << endl; return false; }
        MatrixXd P_eigs = eigenvectors.inverse() * P.transpose(); // Rotate points from world to the coordinate system of the eigenvectors // Pcopy is transposed (i.e. 2xnP)
    	VectorXd extents(eigenvalues.size());
        for (int i = 0; i < eigenvalues.size(); ++i)
        {
            extents[i] = P_eigs.row(i).maxCoeff() - P_eigs.row(i).minCoeff();
        }
        eigenvalues = extents;
    }*/
	// std::cout << "The eigenvalues of P are:\n" << eigenvalues << std::endl;
	// std::cout << "The eigenvectors of P are:\n" << eigenvectors << std::endl;

    return true;
}


