#include "Membrane.h"
#include "General_functions.hpp"

double Membrane::Average_Node_Distance()
{
    double average_membrane_Node_distance=0.0;
    double temp[3];
    double length;
    int Node_A, Node_B;
    for (int i=0; i<Num_of_Triangle_Pairs; i++)
    {
        Node_A=Node_Bond_list[i][0];
        Node_B=Node_Bond_list[i][1];
        temp[0]=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        temp[1]=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        temp[2]=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        length= vector_length(temp);
        average_membrane_Node_distance+=length;
    }
    average_membrane_Node_distance=average_membrane_Node_distance/Num_of_Triangle_Pairs;
    return(average_membrane_Node_distance);
    
}
 
void Membrane::calculate_volume_and_surface_area(){
    volume = 0;
    surface_area=0;
    
    update_COM_position();
    
    int node_A, node_B, node_C;
    double AB[3], AC[3], height[3], ABxAC[3];
    
    for(  int i=0;i<Triangle_list.size();i++)
    {
        node_A = Triangle_list[i][0];
        node_B = Triangle_list[i][1];
        node_C = Triangle_list[i][2];
        
        AB[0] = Node_Position[node_B][0] - Node_Position[node_A][0];
        AB[1] = Node_Position[node_B][1] - Node_Position[node_A][1];
        AB[2] = Node_Position[node_B][2] - Node_Position[node_A][2];
        
        AC[0] = Node_Position[node_C][0] - Node_Position[node_A][0];
        AC[1] = Node_Position[node_C][1] - Node_Position[node_A][1];
        AC[2] = Node_Position[node_C][2] - Node_Position[node_A][2];
        
        height[0] = COM_position[0] - Node_Position[node_A][0];
        height[1] = COM_position[1] - Node_Position[node_A][1];
        height[2] = COM_position[2] - Node_Position[node_A][2];
        
        crossvector(ABxAC, AB, AC);
        
        double area = 0.5*vector_length(ABxAC);
        double H = abs(innerproduct(ABxAC, height)/(2*area) );
        
        surface_area += area;
        volume += H*area/3.0;
        
    }
    
}

double Membrane::calc_theta_angle_ABC(int node_A, int node_B, int node_C){
    /**A is the is the middle point of the angle*/
    double AB[3], AC[3], height[3], ACxAB[3];
    AB[0] = Node_Position[node_B][0] - Node_Position[node_A][0];
    AB[1] = Node_Position[node_B][1] - Node_Position[node_A][1];
    AB[2] = Node_Position[node_B][2] - Node_Position[node_A][2];
   
    AC[0] = Node_Position[node_C][0] - Node_Position[node_A][0];
    AC[1] = Node_Position[node_C][1] - Node_Position[node_A][1];
    AC[2] = Node_Position[node_C][2] - Node_Position[node_A][2];
   
    crossvector(ACxAB, AB, AC);
    return innerproduct(AB, AC)/vector_length(ACxAB);
}

void Membrane::calculate_surface_area_with_voronoi(){
    surface_area=0;
    
    update_COM_position();
    
    int num_5 = 0;
    int num_6 = 0;
    int num = 0;
    for (int i=0; i<Num_of_Nodes; i++){
        int size = Node_neighbour_list[i].size();
        if (size == 6) {
            num_6++;
        } else if (size == 5){
            num_5++;
        } else{
            num++;
        }
    }
    
    cout<<"\nvoronoi check:"<<endl;
    cout<<"nodes with coordiante number 5: "<<num_5<<endl;
    cout<<"nodes with coordiante number 6: "<<num_6<<endl;
    cout<<"nodes with coordiante number other than 5 and 6: "<<num<<endl;
    calculate_volume_and_surface_area();
    cout<<"The surface area calculated using the triangle surface sum: "<<return_surface_area()<<endl;
    
    int triangle1_node_A,
        triangle1_node_B,
        triangle1_node_C,
        triangle2_node_D;
    
    vector<vector<double> > cot_theta_list;
    cot_theta_list.resize(Num_of_Node_Pairs);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        cot_theta_list[i].resize(2,0);
        
        triangle1_node_C = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][0];
        triangle1_node_B = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][1];
        triangle1_node_A = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][2];
        triangle2_node_D = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][3];
        
        cot_theta_list[i][0] = calc_theta_angle_ABC(triangle1_node_C, triangle1_node_A, triangle1_node_B );
        cot_theta_list[i][1] = calc_theta_angle_ABC(triangle2_node_D, triangle1_node_A, triangle1_node_B );
    }
    
    
    cout<<"here\n";
    double sigma_sum=0;

    for (int i=0; i<Num_of_Nodes; i++) {
        node_voronoi_area[i]=0;
        for (int j=0; j<Node_neighbour_list[i].size(); j++) {
//            cout<<"Node_neighbour_list["<<i<<"]["<<j<<"]="<<node_pair_vec[i][j]<<endl;
            int node_1 = i;
            int node_2 = Node_neighbour_list[node_1][j];
            int bond12 = Node_neighbour_list_respective_bond_index[node_1][j];
            
            double bond_vec[3]={Node_Position[node_1][0]-Node_Position[ node_2 ][0],
                                Node_Position[node_1][1]-Node_Position[ node_2 ][1],
                                Node_Position[node_1][2]-Node_Position[ node_2 ][2],
                                };
            node_voronoi_area[i] += 0.125 * vector_length(bond_vec)* vector_length(bond_vec)* (cot_theta_list[bond12 ][0] +
                 cot_theta_list[bond12 ][1]);
        }
        sigma_sum += node_voronoi_area[i];
    }
    cout<<"The surface area calculated using the voronoi surface sum: "<<sigma_sum<<endl;
    cout<<"surface_area-sigmaij = "<<surface_area-sigma_sum<<endl;
    cout<<"\nDone\n";
}
