//
//  ECM_read_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"
#include <sstream>

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

//void ECM::read_gmesh_file_3D(std::string gmesh_file){
//    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
//    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
//    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
//    std::string temp_string;
//    for (int i=0; i<6; i++)
//    {
//        read>> temp_string;
//    }
//    read>> Num_of_Nodes;
//
//    Node_Velocity.resize(Num_of_Nodes);
//    Node_Force.resize(Num_of_Nodes);
//    for(int i=0;i<Num_of_Nodes;i++)
//    {
//        Node_Velocity[i].resize(3,0);
//        Node_Force[i].resize(3,0);
//    }
//
//    // In this section the Node coordinates are read from the Gmesh membrane generated file. These include both the Nodes on the Membrane and on the nucleus membrane.
//    vector<double> temp_node_position;
//    temp_node_position.resize(3);
//    for(int i=0;i<Num_of_Nodes;i++)
//    {
//        read>> temp_int;
//        read>> temp_node_position[0];
//        read>> temp_node_position[1];
//        read>> temp_node_position[2];
//        temp_node_position[0]*=rescale_factor;
//        temp_node_position[1]*=rescale_factor;
//        temp_node_position[2]*=rescale_factor;
//        Node_Position.push_back(temp_node_position);
//    }
//    read>> temp_string;
//    read>> temp_string;
//    read>> temp_int;
//    Num_of_Triangles = temp_int;
//
//    vector<int> push;
//    push.resize(3);
//    for(int i=0;i<Num_of_Triangles;i++)
//    {
//        read>>temp_int;
//        read>>temp_int;
//        read>>temp_int;
//        read>>temp_int;
//        read>>temp_int;
//
//        read>>push[0];
//        read>>push[1];
//        read>>push[2];
//        //We have written the programme so that the Node indecies start from '0'. The node indecies in the Gmesh generated file start from '1', so we correct the 'Membrane_triangle_list'.
//        push[0]--;
//        push[1]--;
//        push[2]--;
//        Triangle_List.push_back(push);
//        // Sometimes Gmesh will create duplicate nodes, this will delete any duplicates
//        if (i!=0)
//        {
//            if (Triangle_List[Triangle_List.size()-1][0]==Triangle_List[Triangle_List.size()-2][0] && Triangle_List[Triangle_List.size()-1][1]==Triangle_List[Triangle_List.size()-2][1] && Triangle_List[Triangle_List.size()-1][2]==Triangle_List[Triangle_List.size()-2][2])
//            {
//                cout<<Triangle_List[Triangle_List.size()-1][0]<<"\t"<<Triangle_List[Triangle_List.size()-1][1]<<"\t"<<Triangle_List[Triangle_List.size()-1][2]<<endl;
//                Triangle_List.erase(Triangle_List.begin()+Triangle_List.size()-2);
//            }
//        }
//    }
//
//}








void ECM::read_gmesh_file_3D(std::string gmesh_file){
    cout<<endl<<endl<<gmesh_file<<endl<<endl;
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    std::string temp_string;

    int Num_of_objects;

    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
    read>> Num_of_Nodes;

    Node_Velocity.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
    Node_Force.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
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

    // In this section the Node list that make up triangles on the outer membrane and nucleus are read from the Gmesh generated file.

    read>> temp_string;
    read>> temp_string;
    read>> Num_of_objects;

    vector<int> push;
    push.resize(4);
    std::string line;

    getline(read, line);//The file reader is now at the end of the line, we need to read an 'empty' line before going to the next line.

    for (int i=0; i<Num_of_objects; i++) {

        getline(read, line);
        std::istringstream iss(line);
        vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());


        if (split.size()==8) {
            continue;
        } else {
            push[0]= stoi(split[5])-1;
            push[1]= stoi(split[6])-1;
            push[2]= stoi(split[7])-1;
            push[3]= stoi(split[8])-1;
            Pyramid_Nodes.push_back(push);
        }
    }
}








void ECM::read_gmesh_file_3D_square(std::string gmesh_file){
    cout<<endl<<endl<<gmesh_file<<endl<<endl;
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    std::string temp_string;

    int Num_of_objects;

    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
    read>> Num_of_Nodes;

    Node_Velocity.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
    Node_Force.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
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

    // In this section the Node list that make up triangles on the outer membrane and nucleus are read from the Gmesh generated file.

    read>> temp_string;
    read>> temp_string;
    read>> Num_of_objects;

    vector<int> push;
    push.resize(8);
    std::string line;

    getline(read, line);//The file reader is now at the end of the line, we need to read an 'empty' line before going to the next line.

    for (int i=0; i<Num_of_objects; i++) {

        getline(read, line);
        std::istringstream iss(line);
        vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());


        if (split.size()==9) {
            continue;
        } else {
            push[0]= stoi(split[5])-1;
            push[1]= stoi(split[6])-1;
            push[2]= stoi(split[7])-1;
            push[3]= stoi(split[8])-1;
            push[4]= stoi(split[9])-1;
            push[5]= stoi(split[10])-1;
            push[6]= stoi(split[11])-1;
            push[7]= stoi(split[12])-1;
            square_Nodes.push_back(push);
        }
    }
}
