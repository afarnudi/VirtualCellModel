#include "Membrane.h"

using std::ifstream;

void Membrane::read_gmesh_file (std::string gmesh_file)
{
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(gmesh_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    std::string temp_string;
    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
    read>> Num_of_Nodes;
    Node_Position.clear();
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
    read>> temp_int;
    Num_of_Triangles=temp_int;
    Triangle_list.clear();
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


void Membrane::import(std::string import_file_name){
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
    shift_position(Shift_position_xyzVector[0], Shift_position_xyzVector[1], Shift_position_xyzVector[2]);
    cout<<"\n\nMembrane class initiated.\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

using namespace std;
void Membrane::read_ply_file (std::string ply_file)
{
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(ply_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary ply file intigers. We never use these intigers in the actual programme.
    std::string temp_string;
    getline(read, temp_string);
    getline(read, temp_string);
    getline(read, temp_string);
    getline(read, temp_string);
    std::istringstream iss(temp_string);
    vector<std::string> words(std::istream_iterator<std::string>{iss}, istream_iterator<std::string>());
    Num_of_Nodes= stoi(words[words.size()-1]);
    
    
    Node_Position.clear();
	Node_Velocity.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
    Node_Force.resize(Num_of_Nodes); //initialize the size of vector witch is just read from "membrane".txt
    for(int i=0;i<Num_of_Nodes;i++)
    {
            Node_Velocity[i].resize(3,0);
            Node_Force[i].resize(3,0);
    }
    getline(read, temp_string);
    getline(read, temp_string);
    getline(read, temp_string);
    getline(read, temp_string);
    std::istringstream iss2(temp_string);
    vector<std::string> words2(std::istream_iterator<std::string>{iss2}, istream_iterator<std::string>());
    Num_of_Triangles =stoi(words2[words2.size()-1]);
    getline(read, temp_string);
    getline(read, temp_string);
    // In this section the Node coordinates are read from the ply file.
    vector<double> temp_node_position;
    temp_node_position.resize(3);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        read>> temp_node_position[0];
        read>> temp_node_position[1];
        read>> temp_node_position[2];
        temp_node_position[0]*=rescale_factor;
        temp_node_position[1]*=rescale_factor;
        temp_node_position[2]*=rescale_factor;
        Node_Position.push_back(temp_node_position);
    }
    
    // In this section the Node list that make up triangles on the outer membrane and nucleus are read from the  ply file.
    Triangle_list.clear();
    vector<int> push;
    push.resize(3);
    for(int i=0;i<Num_of_Triangles;i++)
    {
        read>>temp_int;        
        read>>push[0];
        read>>push[1];
        read>>push[2];
        Triangle_list.push_back(push);
        
    }
    
    
}
