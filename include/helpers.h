//
// Created by Ioanna Mitropoulou on 29.11.21.
//

#pragma once
#ifndef NP_3DP_HELPERS_H
#define NP_3DP_HELPERS_H

#include <igl/matlab_format.h>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <deque>
#include <set>
#include <unsupported/Eigen/SparseExtra>

// Platform-specific includes
#if defined(_WIN32)
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

using namespace Eigen;
using namespace std;

namespace Helpers {

    template <class T>
    inline bool in_vector(const std::vector<T>& vec, const T& item) {
        return std::count(vec.begin(), vec.end(), item);
    }

    template <class T>
    inline bool in_deque(const std::deque<T>& Q, const T& item) {
        return std::find(Q.begin(), Q.end(), item) != Q.end();
    }

    template <class T>
    inline bool in_set(const set<T>& Q, const T& item) {
        return std::find(Q.begin(), Q.end(), item) != Q.end();
    }

    template <class T>
    inline set<T> common_elements_in_sets(const set<T>& s1, const set<T>& s2)
    {
        set<int> intersection;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(intersection, intersection.begin()));
        return intersection;
    }

    template<class T>
    inline vector<T> common_elements_in_vectors(const vector<T>& v1, const vector<T>& v2)
    {
        //vector<T> v3(v1.size() + v2.size());
        //vector<T>::iterator it = set_intersection(v1.begin(), v1.end(),v2.begin(), v2.end(), v3.begin());
        //return v3;

        vector<T> v3;
        for (T vi : v1)
            if (in_vector(v2, vi))
                v3.push_back(vi);
        return v3;
    }

    template <class SomeVector>
    inline bool write_vector_to_txt(const SomeVector& s, const std::string& fileName, bool quiet = false) {
        std::ofstream f(fileName);
        for (int i = 0; i < s.size(); i++)
        {
            f << s[i] << " , ";
        }
        f.close();
        if (!quiet)
            cout << "Saved : " << fileName << endl;
        return !f.fail();
    }

	template <class SomeVector>
    inline bool write_numbered_vector_to_txt(const SomeVector& s, const std::string& fileName, bool quiet = false) {
        std::ofstream f(fileName);
        for (int i = 0; i < s.size(); i++)
        {
            f << i << " , " << s[i] << endl;
        }
        f.close();
        if (!quiet)
            cout << "Saved : " << fileName << endl;
        return !f.fail();
    }

    template<class T>
    inline int element_index(const vector<T>& vec, const T& element) {
        if (!in_vector(vec, element)) { return -1; };
        return std::find(vec.begin(), vec.end(), element) - vec.begin();
    }

    template <class EigenMatrix>
    void removeRow(EigenMatrix& matrix, unsigned int rowToRemove) // removes a row from an Eigen::MatrixXd
    {
        unsigned int numRows = matrix.rows() - 1;
        unsigned int numCols = matrix.cols();

        if (rowToRemove < numRows)
            matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

        matrix.conservativeResize(numRows, numCols);
    }

    template <class EigenMatrix>
    void removeRowPtr(EigenMatrix* matrix, unsigned int rowToRemove) // removes a row from an Eigen::MatrixXd
    {
        unsigned int numRows = matrix->rows() - 1;
        unsigned int numCols = matrix->cols();

        if (rowToRemove < numRows)
            matrix->block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix->block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

        matrix->conservativeResize(numRows, numCols);
    }

    template <class EigenMatrix>
    void removeCol(EigenMatrix& matrix, unsigned int colToRemove) // removes a row from an Eigen::MatrixXd
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols() - 1;

        if (colToRemove < numCols)
            matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

        matrix.conservativeResize(numRows, numCols);
    }

    template <class T>
    inline void print_vector(const std::vector<T>& vec) {
        for (auto& element : vec) { cout << element << endl; }
    }

    template <class SomeEigenMatrix>
    inline bool write_matrix_to_txt(const SomeEigenMatrix& s, const std::string& fileName, bool quiet = false)
    {
        std::ofstream f(fileName);
        for (int i = 0; i < s.rows(); i++)
        {
            for (int j = 0; j < s.cols(); ++j)
            {
                f << s(i, j) << " ";
            }
            f << endl;
        }
        f.close();
        if (!quiet)
            cout << "Saved : " << fileName << endl;
        return !f.fail();
    }

    template <class T>
    inline set<T> convert_vector_to_set(const vector<T>& v) {
        set<T> s;
        for (int x : v)
            s.insert(x);
        return s;
    }

    MatrixXd read_matrix_Xx3_from_txt(const std::string& fileName);
    MatrixXd read_matrix_Xx6_from_txt(const std::string& fileName);

    VectorXi read_vectorXi_from_txt(const std::string& fileName);
    VectorXi read_numbered_vectorXi_from_txt(const std::string& fileName);

    bool read_vf_constraints(const std::string& fileName, VectorXi& fis, MatrixXd& dirs, VectorXd& ws);

    typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> K;
    void setup_line_texture(K& texture_R, K& texture_G, K& texture_B, int texture_scale, int lineWidth);

    double remap_unbound(double input_val, double in_from, double in_to, double out_from, double out_to);
    double remap(double input_val, double in_from, double in_to, double out_from, double out_to);

    inline std::string get_filename(const std::string& folder, const std::string& name, int ID, const std::string& extension) {
        return folder + name + "_id" + std::to_string(ID) + "." + extension;
    }

    template <class T>
    inline void write_to_matlab_format(const T& content, const string& output_folder, const string& name) {
        char filename[512];
        string folder_and_name = output_folder + "/" + name + ".m";
        sprintf(filename, (output_folder + "/" + name + ".m").c_str());
        ofstream ofs;
        ofs.open(filename, std::ofstream::out);
        ofs << igl::matlab_format(content, name);
        ofs.close();
    }

    void write_VV_to_txt(const std::string& fileName, const vector<vector<int>>& VV);

    void write_Emap_to_txt(const string& fileName, const VectorXi& Emap, const MatrixXi& EVquad, bool quiet = false);

    void remove_directory_contents(const std::string& directory_path);
    void create_directory_if_not_exists(const std::string& path);


    // --- deleted igl::diag() functions (last version where they existed was 2.4.0) --->
	template <typename T>
	inline void diag(const Eigen::SparseMatrix<T>& X, Eigen::SparseVector<T>& V){ V = X.diagonal().sparseView();}

	template <typename T,typename DerivedV>
	inline void diag(const Eigen::SparseMatrix<T>& X, Eigen::MatrixBase<DerivedV>& V){ V = X.diagonal(); }

	template <typename T>
	inline void diag(const Eigen::SparseVector<T>& V, Eigen::SparseMatrix<T>& X) {
		// clear and resize output
		Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(V.size(),V.size());
		dyn_X.reserve(V.size());
		// loop over non-zeros
		for(typename Eigen::SparseVector<T>::InnerIterator it(V); it; ++it){
			dyn_X.coeffRef(it.index(),it.index()) += it.value();
		}
		X = Eigen::SparseMatrix<T>(dyn_X);
	}

	template <typename T, typename DerivedV>
	inline void diag(const Eigen::MatrixBase<DerivedV>& V, Eigen::SparseMatrix<T>& X){
		assert(V.rows() == 1 || V.cols() == 1);
		// clear and resize output
		Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(V.size(),V.size());
		dyn_X.reserve(V.size());
		// loop over non-zeros
		for(int i = 0;i<V.size();i++){
			dyn_X.coeffRef(i,i) += V[i];
		}
		X = Eigen::SparseMatrix<T>(dyn_X);
	}
    // <--- igl::diag() functions
}


namespace Geometry
{
    double get_min_distance_from_pts(const Vector3d& pt, const vector<Vector3d>& cloud);
    
    Vector3d project_pt_on_line(const Vector3d& A, const Vector3d& B, const Vector3d& P);
    bool projected_point_inside_segment(const Vector3d& p, const Vector3d& q, const Vector3d& r); // returns true if q is on the segment p-r.
    Vector3d project_pt_on_plane(const Vector3d& N, const Vector3d& origin, const Vector3d& pt);
    pair<Vector3d, bool> intersection_line_line(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d, double tol);
    pair<Vector3d, bool> intersection_line_plane(const Vector3d& a, const Vector3d& b, const Vector3d& o, const Vector3d& n, double tol);
	pair<Vector3d, bool> projected_intersection_line_line(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d, const Vector3d& plane_N, const Vector3d& plane_origin);
    Vector3d project_pt_on_line_segment(const Vector3d& A, const Vector3d& B, const Vector3d& P);
    Matrix4d rotation_matrix_from_vectors(const Vector3d& vec1, const Vector3d& vec2);
    double distance_of_pt_from_plane(const RowVector3d& vertex, const RowVector3d& N, const RowVector3d& C);

    bool PCA(const MatrixXd& P, VectorXd& eigenvalues, MatrixXd& eigenvectors);

}

#endif //NP_3DP_HELPERS_H