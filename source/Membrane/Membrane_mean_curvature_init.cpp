#include "Membrane.h"
#include "General_functions.hpp"

void Membrane::mean_curvature_init(void)
{
    /**
     Create a list for each node in the mesh. The first element is the node index. The rest of the elements contain the indices of nodes that are connected to that node (element zero) with a bond.
     The list of nodes are in order meaning that all consecutive pairs make a triangle with the first node (element zero).
     Example:
     
     List = [0 , 1, 4, 6, 9 , 8]
          0 is the central node and [1, 4, 6, 9, 8] all share a bond with node 0.
          [0, 1, 4], [0, 4, 6], [0, 6, 9], , [0, 9, 8], , [0, 8, 1] make triangles.
     */
    vector<vector<int> > Ordered_Node_neighbour_list;
    Ordered_Node_neighbour_list.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        
        vector<int> neighbour_list = Node_neighbour_list[i];
//        cout<<neighbour_list.size()<<endl;
        Ordered_Node_neighbour_list[i].push_back(i);
        Ordered_Node_neighbour_list[i].push_back(neighbour_list[0]);
        neighbour_list.erase(neighbour_list.begin() );
        while (neighbour_list.size()!=0) {
            int last_node = Ordered_Node_neighbour_list[i][Ordered_Node_neighbour_list[i].size()-1];
//            cout<<"last_node "<<last_node<<"; neighbour_list.size() "<<neighbour_list.size()<<endl;
            for (int j=0; neighbour_list.size(); j++) {
                int node_neighbour = neighbour_list[j];
                
                if (check_if_nodes_make_triangle(i,node_neighbour,last_node)) {
                    Ordered_Node_neighbour_list[i].push_back(node_neighbour);
                    neighbour_list.erase(neighbour_list.begin() + j);
                    break;
                }
                
            }
        }
        
    }
    int max_node_order=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        int node_order = Ordered_Node_neighbour_list[i].size()-1;
        if (max_node_order<node_order) {
            max_node_order = node_order;
        }
    }
//    cout<<max_node_order<<endl;
    nodeOrder_NodeIndex_NodeNeighbourList.resize(max_node_order+1);
    for (int i=0; i<Num_of_Nodes; i++) {
        int node_order = Ordered_Node_neighbour_list[i].size()-1;
        nodeOrder_NodeIndex_NodeNeighbourList[node_order].push_back(Ordered_Node_neighbour_list[i]);
    }
}



bool Membrane::check_if_nodes_make_triangle(int node1,int node2, int node3){
    bool triangle_exists = false;
    
    for (int i=0; i<Triangle_list.size(); i++) {
        
        int tri_node1 = Triangle_list[i][0];
        int tri_node2 = Triangle_list[i][1];
        int tri_node3 = Triangle_list[i][2];
        
        
        if (tri_node1 == node1 || tri_node1 == node2 || tri_node1 == node3) {
            if (tri_node2 == node1 || tri_node2 == node2 || tri_node2 == node3){
                if (tri_node3 == node1 || tri_node3 == node2 || tri_node3 == node3){
                    triangle_exists=true;
//                    cout<<tri_node1<<" "<<tri_node2<<" "<<tri_node3<<" == "<<node1<<" "<<node2<<" "<<node3<<endl;
                    break;
                }
                
            }
        }
    }
    
    
    return triangle_exists;
}
