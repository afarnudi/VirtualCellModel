#include "Chromatin.h"

using std::cout;
using std::endl;

void Chromatin::import_coordinates(string import_file_name){
    cout<<TPINK<<"Importing Chromatin coordinates:"<<endl;
//    cout<<import_file_name<<endl<<endl;
    std::ifstream read_coordinate_file;
    
    read_coordinate_file.open(import_file_name.c_str());
    if ( read_coordinate_file.is_open() ) {
        cout << "\nCoordinate file opened "<<TSUCCESS<<"successfully. \n\n"<<TRESET;
    }else{
        cout << TFAILED<<"Unable to read"<<TRESET<<" '"<<TFILE<<import_file_name<<TRESET<<"' coordinate file.\n";
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
    
    cout<<"Number of nodes\t"<<Num_of_Nodes<<endl;
    
    bond_length/=Num_of_Nodes-1.0;
    cout<<"Average bond length\t"<<bond_length<<endl;
    cout<<TWARN<<"Node forces set to zero"<<TRESET<<endl;
    
    
    
    
    if (rescale_factor!=1 && rescale_factor!=0) {
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]*=rescale_factor;
            Node_Position[i][1]*=rescale_factor;
            Node_Position[i][2]*=rescale_factor;
        }
        bond_length *= rescale_factor;
        cout<<"Average bond length rescaled to: \t"<<bond_length<<endl;
        
    }
    
    if (bond_radius> 0.0001) {
        if (bond_length<2) {
            cout<<"Bond length is shorter than the node diameter. The relative parameters, bond_length and node_radius, can be adjusted in the chromatin configuration file.\n";
            exit(EXIT_FAILURE);
        } else {
            cout<<"will generate virtual spheres in between node beads to prevent bonds from slipping through each other.\n";
            generate_virtual_sites();
        }
    }
    
    ABC_index.resize(Num_of_Nodes);
    pdb_label_check();
    shift_node_positions();
    cout<<TSUCCESS<<"\n\nChromatin class initiated.\n"<<TRESET;
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
    cout<<"Node distances:"<<endl;
    cout<<"Max "<<Max_node_pair_length<<"\tMin "<<Min_node_pair_length<<"\tAverage "<<Average_node_pair_length<<endl;
    
    shift_node_positions();
    pdb_label_check();
    
    cout<<TSUCCESS<<"\n\nChromatin class initiated.\n"<<TRESET;
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}
