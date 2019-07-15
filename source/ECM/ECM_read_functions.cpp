//
//  ECM_read_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"

void ECM::read_gmesh_file_2D(std::string gmesh_file){
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    std::string temp_string;
    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
    read>> Num_of_Nodes;
    
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
    }
    
    // In this section the Node coordinates are read from the Gmesh membrane generated file. These include both the Nodes on the Membrane and on the nucleus membrane.
    vector<double> temp_node_position;
    temp_node_position.resize(3);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        read>> temp_int;
        read>> temp_node_position[0];
        read>> temp_node_position[1];
        read>> temp_node_position[2];
        temp_node_position[0]*=rescale_factor;
        temp_node_position[1]*=rescale_factor;
        temp_node_position[2]*=rescale_factor;
        Node_Position.push_back(temp_node_position);
    }
    read>> temp_string;
    read>> temp_string;
    read>> temp_int;
    Num_of_Triangles = temp_int;
    
    vector<int> push;
    push.resize(3);
    for(int i=0;i<Num_of_Triangles;i++)
    {
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        
        read>>push[0];
        read>>push[1];
        read>>push[2];
        //We have written the programme so that the Node indecies start from '0'. The node indecies in the Gmesh generated file start from '1', so we correct the 'Membrane_triangle_list'.
        push[0]--;
        push[1]--;
        push[2]--;
        Triangle_List.push_back(push);
        // Sometimes Gmesh will create duplicate nodes, this will delete any duplicates
        if (i!=0)
        {
            if (Triangle_List[Triangle_List.size()-1][0]==Triangle_List[Triangle_List.size()-2][0] && Triangle_List[Triangle_List.size()-1][1]==Triangle_List[Triangle_List.size()-2][1] && Triangle_List[Triangle_List.size()-1][2]==Triangle_List[Triangle_List.size()-2][2])
            {
                cout<<Triangle_List[Triangle_List.size()-1][0]<<"\t"<<Triangle_List[Triangle_List.size()-1][1]<<"\t"<<Triangle_List[Triangle_List.size()-1][2]<<endl;
                Triangle_List.erase(Triangle_List.begin()+Triangle_List.size()-2);
            }
        }
    }
    
}

void ECM::read_gmesh_file_3D(std::string gmesh_file){
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    std::string temp_string;
    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
    read>> Num_of_Nodes;
    
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
    }
    
    // In this section the Node coordinates are read from the Gmesh membrane generated file. These include both the Nodes on the Membrane and on the nucleus membrane.
    vector<double> temp_node_position;
    temp_node_position.resize(3);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        read>> temp_int;
        read>> temp_node_position[0];
        read>> temp_node_position[1];
        read>> temp_node_position[2];
        temp_node_position[0]*=rescale_factor;
        temp_node_position[1]*=rescale_factor;
        temp_node_position[2]*=rescale_factor;
        Node_Position.push_back(temp_node_position);
    }
    read>> temp_string;
    read>> temp_string;
    read>> temp_int;
    Num_of_Triangles = temp_int;
    
    vector<int> push;
    push.resize(3);
    for(int i=0;i<Num_of_Triangles;i++)
    {
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        
        read>>push[0];
        read>>push[1];
        read>>push[2];
        //We have written the programme so that the Node indecies start from '0'. The node indecies in the Gmesh generated file start from '1', so we correct the 'Membrane_triangle_list'.
        push[0]--;
        push[1]--;
        push[2]--;
        Triangle_List.push_back(push);
        // Sometimes Gmesh will create duplicate nodes, this will delete any duplicates
        if (i!=0)
        {
            if (Triangle_List[Triangle_List.size()-1][0]==Triangle_List[Triangle_List.size()-2][0] && Triangle_List[Triangle_List.size()-1][1]==Triangle_List[Triangle_List.size()-2][1] && Triangle_List[Triangle_List.size()-1][2]==Triangle_List[Triangle_List.size()-2][2])
            {
                cout<<Triangle_List[Triangle_List.size()-1][0]<<"\t"<<Triangle_List[Triangle_List.size()-1][1]<<"\t"<<Triangle_List[Triangle_List.size()-1][2]<<endl;
                Triangle_List.erase(Triangle_List.begin()+Triangle_List.size()-2);
            }
        }
    }
    
}
