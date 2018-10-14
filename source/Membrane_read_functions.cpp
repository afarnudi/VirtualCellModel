#include "Membrane.h"

void Membrane::read_gmesh_file (string gmesh_file)
{
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;
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
    read>> temp_int;
    Num_of_Triangles=temp_int;
    
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
        //We have written the programme so that the Node indecies start from '0'. The node indecies in the Gmesh generated file start from '1', so we correct the 'Triangle_list'.
        push[0]--;
        push[1]--;
        push[2]--;
        Triangle_list.push_back(push);
        // Sometimes Gmesh will create duplicate nodes, this will delete any duplicates
        if (i!=0)
        {
            if (Triangle_list[Triangle_list.size()-1][0]==Triangle_list[Triangle_list.size()-2][0] && Triangle_list[Triangle_list.size()-1][1]==Triangle_list[Triangle_list.size()-2][1] && Triangle_list[Triangle_list.size()-1][2]==Triangle_list[Triangle_list.size()-2][2])
            {
                cout<<Triangle_list[Triangle_list.size()-1][0]<<"\t"<<Triangle_list[Triangle_list.size()-1][1]<<"\t"<<Triangle_list[Triangle_list.size()-1][2]<<endl;
                Triangle_list.erase(Triangle_list.begin()+Triangle_list.size()-2);
            }
        }
    }
}


void Membrane::read_membrabe_input(string input_file)

{
    /*start of initializing constants*/
    ifstream inputs;
    inputs.open("membrane_inputs_file_name");
    string temp_str; //This is just a temp string charachter that we use to read unnecessary words in inputs file . We never use this  in the actual programme.
    inputs>>temp_str;
    inputs>>Spring_coefficient; // streching constant
    inputs>>temp_str;
    inputs>>Bending_coefficient; // bending constant
    inputs>>temp_str;
    inputs>>Damping_coefficient; // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
    inputs>>temp_str;
    inputs>>Node_Mass;
    inputs>>temp_str;
    inputs>>K_surfaceConstant_local;
    inputs>>temp_str;
    inputs>>on_or_off_Spring_force_cutt_off; //??? I add it myself because virus should not have cut off
    inputs>>temp_str;
    inputs>>Spring_force_cutt_off;
    inputs>>temp_str;
    inputs>>ShiftinXdirection; //???
    inputs>>temp_str;
    inputs>>ShiftinZdirection; //???
    inputs>>temp_str;
    inputs>>Downward_speed; //???
    inputs>>temp_str;
    //inputs>>energy_calculation_flag;
    /*end of initializing constants*/
}

void Membrane::import(string import_file_name){
    cout<<"Importing the Membrane from the resume file:"<<endl;
    cout<<import_file_name<<endl<<endl;
    ifstream read_resume_file;
    
    read_resume_file.open(import_file_name.c_str());
    if (read_resume_file) {
        cout << "Managed to read file successfully. \n\n";
    }else{
        cout << "Unable to read file.";
    }
    int temp;
    read_resume_file>>temp;
    cout<<"Resuming the MD from step "<<temp<<endl;
    
    read_resume_file>>Num_of_Nodes;
    cout<<"Number of Nodes: "<<Num_of_Nodes<<endl;
    
    Node_Force.resize(Num_of_Nodes);
    
    vector<double> read_double;
    read_double.resize(3);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        read_resume_file>>read_double[0]>>read_double[1]>>read_double[2];
        Node_Position.push_back(read_double);

        read_resume_file>>read_double[0]>>read_double[1]>>read_double[2];
        Node_Velocity.push_back(read_double);
        
        Node_Force[i].resize(3,0);
    }
    cout<<"Coordinates and velocities loaded"<<endl;
    cout<<"Node forces set to zero"<<endl;
    
    read_resume_file>>Num_of_Triangles;
    cout<<"Number of Triangles: "<<Num_of_Triangles<<endl;
    
    vector<int> read_int;
    read_int.resize(3);
    for (int i=0; i<Num_of_Triangles; i++) {
        read_resume_file>>read_int[0]>>read_int[1]>>read_int[2];
        Triangle_list.push_back(read_int);
    }
    
    read_resume_file>>Num_of_Node_Pairs;
    cout<<"Number of Node Pairs: "<<Num_of_Node_Pairs<<endl;
    
    
    read_int.resize(2);
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        read_resume_file>>read_int[0]>>read_int[1];
        Node_Bond_list.push_back(read_int);
    }
    //In the import function we should call the neighbour list constructor
    read_resume_file>>Num_of_Triangle_Pairs;
    cout<<"Number of Triangle Pairs: "<<Num_of_Triangle_Pairs<<endl;
    
    vector<int> read_int_4;
    read_int_4.resize(4);
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        read_resume_file>>read_int[0]>>read_int[1];
        Triangle_pair_list.push_back(read_int);
        
        read_resume_file>>read_int_4[0]>>read_int_4[1]>>read_int_4[2]>>read_int_4[3];
        Triangle_Pair_Nodes.push_back(read_int_4);
    }
    
    Node_neighbour_list_constructor();
    read_resume_file>>Max_node_pair_length>>Min_node_pair_length>>Average_node_pair_length;
    cout<<"Spring coefficients:"<<endl;
    cout<<"Max="<<Max_node_pair_length<<"\tmin="<<Min_node_pair_length<<"\tAverage="<<Average_node_pair_length<<endl;
    cout<<"\n\nMembrane class initiated.\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}
