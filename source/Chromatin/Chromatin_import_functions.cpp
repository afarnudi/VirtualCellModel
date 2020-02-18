#include "Chromatin.h"

using std::cout;
using std::endl;

void Chromatin::import_coordinates(string import_file_name){
//    cout<<"Importing the Chromatin from the resume file:"<<endl;
//    cout<<import_file_name<<endl<<endl;
    std::ifstream read_coordinate_file;
    
    read_coordinate_file.open(import_file_name.c_str());
    if ( read_coordinate_file.is_open() ) {
        cout << "Coordinate file opened successfully successfully. \n\n";
    }else{
        cout << "Unable to read coordinate file.\n";
        exit(EXIT_FAILURE);
    }
    
    bond_length=0;
    Num_of_Nodes=0;
    vector<double> read_double;
    read_double.resize(3);
    
    vector<double> zero_double;
    zero_double.resize(3,0);
    
    double temp_vector[3];
    
    while (true) {
        read_coordinate_file>>read_double[0]>>read_double[1]>>read_double[2];
        if (read_coordinate_file.eof()) break;
        Node_Position.push_back(read_double);
        if (Num_of_Nodes>0) {
            temp_vector[0] = Node_Position[Num_of_Nodes][0]-Node_Position[Num_of_Nodes-1][0];
            temp_vector[1] = Node_Position[Num_of_Nodes][1]-Node_Position[Num_of_Nodes-1][1];
            temp_vector[2] = Node_Position[Num_of_Nodes][2]-Node_Position[Num_of_Nodes-1][2];
            
            bond_length += vector_length(temp_vector);
        }
        
        read_coordinate_file>>read_double[0]>>read_double[1]>>read_double[2];
        Node_Velocity.push_back(read_double);
        Node_Force.push_back(zero_double);
        
        Num_of_Nodes++;
    }
    
    cout<<"# of nodes\t"<<Num_of_Nodes<<endl;
    
    bond_length/=Num_of_Nodes-1.0;
    cout<<"bond length set to \t"<<bond_length<<endl;
    
    
    ABC_index.resize(Num_of_Nodes);
    cout<<"Coordinates and velocities loaded"<<endl;
    cout<<"Node forces set to zero"<<endl;
    
    shift_node_positions();
    pdb_label_check();
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Position[i][0]*=rescale_factor;
        Node_Position[i][1]*=rescale_factor;
        Node_Position[i][2]*=rescale_factor;
    }
    
    cout<<"\n\nChromatin class initiated.\n";
}

void Chromatin::import_resume(string import_file_name){
//    cout<<"Importing the Chromatin from the resume file:"<<endl;
//    cout<<import_file_name<<endl<<endl;
    std::ifstream read_resume_file;
    
    read_resume_file.open(import_file_name.c_str());
    if ( read_resume_file.is_open() ) {
        cout << "Resume file opened successfully successfully. \n\n";
    }else{
        cout << "Unable to read file.\n";
        exit(EXIT_FAILURE);
    }
    int temp;
    read_resume_file>>temp;
    cout<<"Resuming the MD from "<<temp<<" Ps"<<endl;
    
    read_resume_file>>Num_of_Nodes;
    cout<<"Number of Nodes: "<<Num_of_Nodes<<endl;
    
    read_resume_file>>num_of_node_types;
    cout<<"Number of Node types: "<<num_of_node_types<<endl;
    
    epsilon_LJ.resize(num_of_node_types);
    sigma_LJ.resize(num_of_node_types);
    for (int i=0; i<num_of_node_types; i++) {
        read_resume_file>>epsilon_LJ[i];
        read_resume_file>>sigma_LJ[i];
    }
    
    
    Node_Force.resize(Num_of_Nodes);
    ABC_index.resize(Num_of_Nodes);
//    Contact_Matrix.resize(Num_of_Nodes);
    
    vector<double> read_double;
    read_double.resize(3);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        
        read_resume_file>>read_double[0]>>read_double[1]>>read_double[2];
        Node_Position.push_back(read_double);

        read_resume_file>>read_double[0]>>read_double[1]>>read_double[2];
        Node_Velocity.push_back(read_double);
        
        read_resume_file>>ABC_index[i];
        
        Node_Force[i].resize(3,0);
//        Contact_Matrix[i].resize(Num_of_Nodes,0);
    }
    cout<<"Coordinates and velocities loaded"<<endl;
    cout<<"Node forces set to zero"<<endl;
    
//    Node_neighbour_list_constructor();
    read_resume_file>>Max_node_pair_length>>Min_node_pair_length>>Average_node_pair_length;
    cout<<"node distance statistics:"<<endl;
    cout<<"Max="<<Max_node_pair_length<<"\tmin="<<Min_node_pair_length<<"\tAverage="<<Average_node_pair_length<<endl;
    
    shift_node_positions();
    pdb_label_check();
    
    cout<<"\n\nMembrane class initiated.\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}
