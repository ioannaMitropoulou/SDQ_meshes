//
//  Created by Ioanna Mitropoulou on 25.11.21.
//

#include "Np3dpContext.hpp"
#include "quad_mesh.h"

#include <directional/curl_matching.h>
#include <directional/polygonal_edge_topology.h>
#include <igl/serialize.h>
#include <igl/writeOBJ.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>



Np3dpContext::Np3dpContext(const std::string& data_folder) : data_folder(data_folder) {
    output_folder = data_folder + string("output/");
    serialize_folder = output_folder + string("serialize/");


    // --- load triangle mesh
    bool found = directional::readOBJ(DATA_PATH + data_folder + "mesh.obj", mesh);
    if (!found) { cerr << " --- Attention! Could not read " << DATA_PATH + data_folder + "mesh.obj" << ". Either the file is missing, or it has some problem." << endl; return; }
    cout << "Loaded triangle mesh : " << DATA_PATH + data_folder + "mesh.obj" << endl;


    // --- read mesh_scale_factor.txt
    {
        ifstream data(DATA_PATH + data_folder + "scale_factor.txt");
        if (data.is_open()) {
            string s;
            while (getline(data, s)) {
                stringstream line(s);
                line >> mesh_scale_factor;
                break;
            }
            data.close();
        }
        else { cerr << " --- Attention! Could not find 'scale_factor.txt' in the data folder : " << DATA_PATH + data_folder << endl; found_correct_files = false; return; }

        double div = 100 * 1.0;// Config::desired_avg_layer_height; // for layer height 1.0 we divide by 100, for h we divide by 100*h
        parametrization_global_scale = mesh_scale_factor / div;
    }


    // --- create output and serialize folders (if they don't already exist)
    {
		// Check and create output folder
        std::string output_folder_path = DATA_PATH + output_folder;
        Helpers::create_directory_if_not_exists(output_folder_path);
        std::cout << "Created folder " << output_folder_path << std::endl;

        std::string serialize_folder_path = DATA_PATH + serialize_folder;
        Helpers::create_directory_if_not_exists(serialize_folder_path);
        std::cout << "Created folder " << serialize_folder_path << std::endl;
    }

    found_correct_files = true;

    quad_mesh = QuadMesh(output_folder);
}


void Np3dpContext::save_all_pieces_data(bool quiet) const
{
    for (int id = 0; id < nP(); ++id) {

        // get edge topology		
        MatrixXi EV, FE, EF, EFi; MatrixXd FEs; VectorXi innerEdges;
        hedra::polygonal_edge_topology(pieces[id].D, pieces[id].F, EV, FE, EF, EFi, FEs, innerEdges);

        // save files
        std::string filename_off = Helpers::get_filename(DATA_PATH + output_folder, "quad_mesh",  id, "off");
        std::string filename_obj = Helpers::get_filename(DATA_PATH + output_folder, "quad_mesh",  id, "obj");
        std::string filename_V = Helpers::get_filename(DATA_PATH + output_folder, "V",  id, "txt");
        std::string filename_F = Helpers::get_filename(DATA_PATH + output_folder, "F",  id, "txt");
        std::string filename_Emap = Helpers::get_filename(DATA_PATH + output_folder, "Emap",  id, "txt");
        std::string filename_Eribs = Helpers::get_filename(DATA_PATH + output_folder, "ERibs",  id, "txt");
        std::string filename_Ftype = Helpers::get_filename(DATA_PATH + output_folder, "Ftype",  id, "txt");

    	hedra::polygonal_write_OFF(filename_off, pieces[id].V, pieces[id].D, pieces[id].F);
        igl::writeOBJ(filename_obj, pieces[id].V, pieces[id].F);
		Helpers::write_matrix_to_txt(pieces[id].V, filename_V, quiet);
		Helpers::write_matrix_to_txt(pieces[id].F, filename_F, quiet);
    	Helpers::write_Emap_to_txt(filename_Emap, pieces[id].Emap, EV, quiet);
        Helpers::write_Emap_to_txt(filename_Eribs, pieces[id].ERibMap, EV, quiet);
        if (!quiet) std::cout << "Saved partitioned meshes : " << filename_off << " , " << filename_obj << endl;
    }
}


void Np3dpContext::invalidate_parametrization_and_quad_related_data()
{
    is_cut_and_combed = false;
    has_uv_coords = false;
    meshCut = directional::TriMesh();
    UVcoords.setZero(0, 0);

    has_quad_mesh = false;
    quad_mesh = QuadMesh(output_folder);
    quad_mesh_is_subdivided = 0;

    Fmap.setZero(0);
    pieces.clear();
    piece_ID = -1;

    PartitioningCuts.setZero(0);
    manually_selected_partition_eis.clear();
}


void Np3dpContext::serialize() const {
	igl::serialize(parametrization_global_scale, DATA_PATH + serialize_folder + "np3dp.vf_scale_factor.igl");
    igl::serialize(is_cut_and_combed, DATA_PATH + serialize_folder + "np3dp.is_cut_and_combed.igl");
    igl::serialize(texture_scale, DATA_PATH + serialize_folder + "np3dp.texture_scale.igl");
    igl::serialize(texture_line_width, DATA_PATH + serialize_folder + "np3dp.texture_line_width.igl");
    igl::serialize(has_uv_coords, DATA_PATH + serialize_folder + "np3dp.has_uv_coords.igl");

    if (has_uv_coords)
    {
        igl::writeOBJ(DATA_PATH + serialize_folder + "np3dp.meshCut.obj", meshCut.V, meshCut.F);
        igl::serialize(UVcoords, DATA_PATH + serialize_folder + "np3dp.UVcoords.igl");
    }

    igl::serialize(has_quad_mesh, DATA_PATH + serialize_folder + "np3dp.has_quad_mesh.igl");
    if (has_quad_mesh)
    {
	    igl::serialize(quad_mesh_is_subdivided, DATA_PATH + serialize_folder + "np3dp.quad_mesh_is_subdivided.igl");
        quad_mesh.serialize(serialize_folder);
    }

    igl::serialize(has_partitioned_quad_mesh, DATA_PATH + serialize_folder + "np3dp.has_partitioned_quad.igl");
    if (has_partitioned_quad_mesh)
    {
	    igl::serialize(Fmap, DATA_PATH + serialize_folder + "np3dp.Fmap.igl");
        for (int pi=0; pi<nP(); ++pi)
        {
            pieces[pi].serialize(pi, DATA_PATH + serialize_folder);
        }
        igl::serialize(PartitioningCuts, DATA_PATH + serialize_folder + "np3dp.PartitioningCuts.igl");
        igl::serialize(PartitioningCuts, DATA_PATH + serialize_folder + "np3dp.PartitioningCuts.igl");
        igl::serialize(manually_selected_partition_eis, DATA_PATH + serialize_folder  + "np3dp.manually_selected_partition_eis_0.igl");
    }
}


void Np3dpContext::deserialize() {
	bool found = directional::readOBJ(DATA_PATH + data_folder + "mesh.obj", mesh);
    if (!found) { cerr << " --- Attention! Could not read " << DATA_PATH + data_folder + "mesh.obj" << ". Either the file is missing, or it has some problem." << endl; return; }
    cout << "Loaded triangle mesh : " << DATA_PATH + data_folder + "mesh.obj" << endl;
    found_correct_files = true;

	igl::deserialize(parametrization_global_scale, DATA_PATH + serialize_folder + "np3dp.vf_scale_factor.igl");
    igl::deserialize(is_cut_and_combed, DATA_PATH + serialize_folder + "np3dp.is_cut_and_combed.igl");
    igl::deserialize(texture_scale, DATA_PATH + serialize_folder + "np3dp.texture_scale.igl");
    igl::deserialize(texture_line_width, DATA_PATH + serialize_folder + "np3dp.texture_line_width.igl");
    igl::deserialize(has_uv_coords, DATA_PATH + serialize_folder + "np3dp.has_uv_coords.igl");

    if (has_uv_coords)
    {
        MatrixXd V; MatrixXi F;
        igl::readOBJ(DATA_PATH + serialize_folder + "np3dp.meshCut.obj", V, F);
        meshCut.set_mesh(V, F);
        igl::deserialize(UVcoords, DATA_PATH + serialize_folder + "np3dp.UVcoords.igl");
    }

    igl::deserialize(has_quad_mesh, DATA_PATH + serialize_folder + "np3dp.has_quad_mesh.igl");
    if (has_quad_mesh)
    {
	    igl::deserialize(quad_mesh_is_subdivided, DATA_PATH + serialize_folder + "np3dp.quad_mesh_is_subdivided.igl");
        quad_mesh.deserialize(serialize_folder);
    }

    igl::deserialize(has_partitioned_quad_mesh, DATA_PATH + serialize_folder + "np3dp.has_partitioned_quad.igl");
    if (has_partitioned_quad_mesh)
    {
	    igl::deserialize(Fmap, DATA_PATH + serialize_folder + "np3dp.Fmap.igl");
        for (int pi=0; pi<nP(); ++pi)
        {
            pieces[pi].deserialize(pi, DATA_PATH + serialize_folder);
        }
        igl::deserialize(PartitioningCuts, DATA_PATH + serialize_folder + "np3dp.PartitioningCuts.igl");
        igl::deserialize(PartitioningCuts, DATA_PATH + serialize_folder + "np3dp.PartitioningCuts.igl");
        igl::deserialize(manually_selected_partition_eis, DATA_PATH + serialize_folder  + "np3dp.manually_selected_partition_eis_0.igl");
    }

}



