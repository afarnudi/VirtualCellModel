#include "Actin.h"
#include <sstream>

using std::cout;
using std::endl;


void Actin::read_gmesh_file (string gmesh_file)
{
//    cout<<endl<<endl<<gmesh_file<<endl<<endl;
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;
    
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
    string line;
    
    getline(read, line);//The file reader is now at the end of the line, we need to read an 'empty' line before going to the next line.
    
    for (int i=0; i<Num_of_objects; i++) {
        
        getline(read, line);
        std::istringstream iss(line);
        vector<string> split(std::istream_iterator<string>{iss}, std::istream_iterator<string>());
        
        
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

//for type 2 actin
void Actin::read_gmesh_file_2 (string gmesh_file)
{
//    cout<<endl<<endl<<gmesh_file<<endl<<endl;
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;

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

    int type;
    vector<int> push;
    push.resize(2);
    string line;

    getline(read, line);//The file reader is now at the end of the line, we need to read an 'empty' line before going to the next line.

    for (int i=0; i<Num_of_objects; i++) {

        getline(read, line);
        std::istringstream iss(line);
        vector<string> split(std::istream_iterator<string>{iss}, std::istream_iterator<string>());
            type = stoi(split[1]);
            push[0]= stoi(split[2])-1;
            push[1]= stoi(split[3])-1;
        
        if(type==1){
            filaments.push_back(push);
        }
        if(type==2){
            abps.push_back(push);
        }
        if(type==3){
            MTs.push_back(push);
        }
    }

    //cout<<filaments.size()<<'\n';
    //cout<<filaments[2][0]<<"and"<<filaments[2][1]<<'\n';
}

void Actin::read_actin_file(void){
    std::ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(Mesh_file_name.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;
    
    read>> temp_string;
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

    read>> temp_string;
    read>> Num_of_Node_Pairs;
    
    vector<int> temp_node_pair;
    temp_node_pair.resize(2);
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        read>> temp_int;
        read>> temp_node_pair[0];
        read>> temp_node_pair[1];
        temp_node_pair[0]--;
        temp_node_pair[1]--;
        
        Node_Bond_list.push_back(temp_node_pair);
    }
    
}
