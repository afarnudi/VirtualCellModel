#include "Chromatin.h"


void Chromatin::import(string import_file_name){
    cout<<"Importing the Chromatin from the resume file:"<<endl;
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
    
    vector<int> read_int;
    read_int.resize(3);
    
    read_int.resize(2);
    
    
//    Node_neighbour_list_constructor();
    read_resume_file>>Max_node_pair_length>>Min_node_pair_length>>Average_node_pair_length;
    cout<<"Spring coefficients:"<<endl;
    cout<<"Max="<<Max_node_pair_length<<"\tmin="<<Min_node_pair_length<<"\tAverage="<<Average_node_pair_length<<endl;
    cout<<"\n\nMembrane class initiated.\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}
