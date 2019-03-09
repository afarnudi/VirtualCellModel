#include "Actin.h"
#include "General_constants.h"
using namespace std;

void Actin::Node_Bond_relaxed_length_initialiser(void){
    double node_pair_distance[3];
    int node_1,node_2;
    
    Node_Bond_relaxed_length.resize(Num_of_Node_Pairs);
    
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        node_1=Node_Bond_list[i][0];
        node_2=Node_Bond_list[i][1];
        
        node_pair_distance[0]=Node_Position[node_2][0]-Node_Position[node_1][0];
        node_pair_distance[1]=Node_Position[node_2][1]-Node_Position[node_1][0];
        node_pair_distance[2]=Node_Position[node_2][2]-Node_Position[node_1][0];
        
        Node_Bond_relaxed_length[i]=vector_length(node_pair_distance);
    }
}
