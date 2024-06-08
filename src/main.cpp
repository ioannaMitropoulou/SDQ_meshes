#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <iostream>
#include "Np3dpContext.hpp"
#include "vector_field.h"
#include "gui.h"

int main(int argc, char* argv[]) {

	// --- unpack input arguments
	if (argc < 2) {
		std::cout << "Incomplete inputs. Provide name of the folder in data as 1st argument." << std::endl;
		return 0;
	}
	const std::string data_folder = std::string(argv[1]);
 
    // --- create mp3dpContext
    shared_ptr<Np3dpContext> np3dp = make_shared<Np3dpContext>(data_folder);
    if (!np3dp->found_correct_files){std::cout << "Could not find all necessary files in folder " << DATA_PATH +  data_folder << std::endl;  return 0;}
 
    // --- create vector field
    shared_ptr<VectorField> vf = make_shared<VectorField>(np3dp);
 
    // --- run gui
    Gui gui = Gui(np3dp, vf);
}
