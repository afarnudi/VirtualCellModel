#include "Actin.h"
#include "General_constants.h"
using namespace std;

void Actin::set_bond_nominal_length(void){
    Node_Bond_Nominal_Length_in_Nm.clear();
    Node_Bond_Nominal_Length_in_Nm.resize(Num_of_Node_Pairs,0);
    if (Node_Bond_Nominal_Length_stat=="Au") {
        
        for (int i=0; i<Num_of_Node_Pairs; i++) {
            Node_Bond_Nominal_Length_in_Nm[i] = sqrt( (Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0]) + (Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1]) + (Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        }
        cout<<"Using mesh initial distances as the springs nominal length."<<endl;
    } else if (Node_Bond_Nominal_Length_stat=="Av"){
        cout<<"Using the average mesh bond distances ("<<Average_node_pair_length<<") as the springs nominal length."<<endl;
        for (int i=0; i<Num_of_Node_Pairs; i++) {
            Node_Bond_Nominal_Length_in_Nm[i] = Average_node_pair_length;
        }
    } else {
        cout<<"Using "<<Node_Bond_user_defined_Nominal_Length_in_Nm<<" as the springs nominal length."<<endl;
        for (int i=0; i<Num_of_Node_Pairs; i++) {
            Node_Bond_Nominal_Length_in_Nm[i] = Node_Bond_user_defined_Nominal_Length_in_Nm;
        }
    }
    
    
}

