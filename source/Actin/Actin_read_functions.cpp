#include "Actin.h"
#include <sstream>

void Actin::read_gmesh_file (string gmesh_file)
{
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
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
        istringstream iss(line);
        vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
        
        
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
