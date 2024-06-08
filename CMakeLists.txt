cmake_minimum_required(VERSION 3.16)
project(np3dp-v2)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Add your project files
file(GLOB SRCFILES
        src/*.cpp
        )
file(GLOB HDRFILES
        include/*.hpp include/*.h
        )

include_directories(
        include/
)


# --------- libigl
include(libigl)
igl_include(glfw)
igl_include(imgui)
igl_include(opengl)
igl_include(restricted matlab)
igl_include(copyleft cgal)

list(APPEND LIBRARIES igl::core igl::opengl igl::imgui  igl::glfw  igl_restricted::matlab igl_copyleft::cgal)


# --------- Directional library
include(Directional)

# --------- libhedra
option(LIBHEDRA_WITH_CERES      "Use Ceres"         OFF)
include(libhedra)
include_directories(${CMAKE_SOURCE_DIR}/external/libhedra/include)


# --------------
include_directories( ${INCLUDE_DIRS} )
add_executable(${PROJECT_NAME} ${SRCFILES} ${HDRFILES})
target_link_libraries(${PROJECT_NAME} ${LIBRARIES})

#set_property(TARGET np3dp-v2 PROPERTY CXX_STANDARD 17) # needs to be after the 'add_exectuable' line, because that's where the TARGET is created
set(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo -framework Accelerate -framework GLUT ")


# ------------- set the value in DATA_PATH 
set(DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data/ CACHE PATH "location of data")
target_compile_definitions(np3dp-v2 PUBLIC "-DDATA_PATH=\"${DATA_PATH}\"")

igl_copy_dll(${PROJECT_NAME})
