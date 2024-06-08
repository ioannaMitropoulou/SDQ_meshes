//
//  Created by Ioanna Mitropoulou
//

#include "structs.h"



//////////////////////////////
// --- AlignedVertices
//////////////////////////////

void AlignedVertices::append_visualization_data(const MatrixXd& Vquad, MatrixXd& P1, MatrixXd& P2, vector<int>& vis_out) const
{
    for (int k=0; k<N; ++k)
    {
        if (aligned[k])
        {
            const vector<int>& path = all_mesh_paths[k];

            int n1 = P1.rows();
            int n2 = P2.rows();
            P1.conservativeResize(n1 + path.size() - 1, 3);
            P2.conservativeResize(n2 + path.size() - 1, 3);
            for (int i = 0; i < path.size(); ++i)
            {
                if (i == 0) {
                    P1.row(n1 + i) = Vquad.row(path[i]);
                }
                else if (i == path.size() - 1) {
                    P2.row(n2 + i - 1) = Vquad.row(path[i]);
                }
                else {
                    P1.row(n1 + i) = Vquad.row(path[i]);
                    P2.row(n2 + i - 1) = Vquad.row(path[i]);
                }
            }

            vis_out.push_back(vis[0]);
            vis_out.push_back(vis[1]);
        }
    }
}




////////////////////////
/// --- StripNetwroks
////////////////////////


const GraphVertexProps& StripNetworks::strip_index_to_properties(int si) const
{
    int strip_dir = StoD[si];
    const StripGraph& G = strip_dir == 0 ? Dgraph : cDgraph;
    int v = StoG[si]; // graph node index
    //return get(boost::vertex_attribute, G)[v];
    return G[v];
}