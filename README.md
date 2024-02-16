# SDQ_meshes

This repository contains the code used for the publication **Fabrication-aware strip-decomposable quadrilateral meshes**, by Ioanna Mitropoulou, Amir Vaxman, Olga Diamanti, Benjamin Dillenburger, published in Computer-Aided Design, Volume 168, 2024, [doi link](https://doi.org/10.1016/j.cad.2023.103666).


In particular, the code implements;
* The vector field optimization method presented in that paper, where two transversal 2-fields are optimized jointly and integrated in two coupled strip networks (publication Section 4, see `vector_field.h`).
* The overlay of the two networks into a Strip-Decomposable Quad (SDQ) mesh (publication Section 4.4, see `quad_mesh.h`).
* Strip-based editing operations to fix topologic defects of the strip networks (publication Section 5, see `quad_mesh_editing.h` and `strips_topology.h`).
* Tracing paths on the strip network (publication Section 6.1, see `paths.h` and `paths_tracing.h`)


*** Description ***
- Input command line arguments
- Where output data is saved
- How to use, basic overview and UI description.
- Grasshopper visualization and txt files. How to use it.
