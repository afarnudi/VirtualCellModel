#include <sstream>

#include "Chromatin.h"

using namespace std;
void Chromatin::import_coordinates(string import_file_name){
    cout<<TPINK<<"Importing Chromatin coordinates:"<<endl;
    
    if (limit_import!=0) {
        cout<<"will only import the first "<<TFILE<<limit_import<<TRESET<<" coordinates."<<endl;
    }
    
    
    std::ifstream read_coordinate_file;
    
    read_coordinate_file.open(import_file_name.c_str());
    if ( read_coordinate_file.is_open() ) {
        cout << "Coordinate file opened "<<TSUCCESS<<"successfully. \n\n"<<TRESET;
    }else{
        cout << TFAILED<<"Unable to read"<<TRESET<<" '"<<TFILE<<import_file_name<<TRESET<<"' coordinate file.\n";
        exit(EXIT_FAILURE);
    }
    
    
    Num_of_Nodes=0;
    vector<double> read_double;
    read_double.resize(3);
    
    vector<double> zero_double;
    zero_double.resize(3,0);
    
    string line;
    while (!read_coordinate_file.eof()) {
        getline(read_coordinate_file, line);
        istringstream iss(line);
        vector<string> words(istream_iterator<string>{iss}, istream_iterator<string>());
        if (words.size()!=6) {
            break;
        }
        read_double[0]=stod(words[0]);
        read_double[1]=stod(words[1]);
        read_double[2]=stod(words[2]);
        
        Node_Position.push_back(read_double);
        
        read_double[0]=stod(words[3]);
        read_double[1]=stod(words[4]);
        read_double[2]=stod(words[5]);
        
        Node_Velocity.push_back(read_double);
        Node_Force.push_back(zero_double);
        
        Num_of_Nodes++;
        if (Num_of_Nodes==limit_import) {
            break;
        }
    }
    cout<<"Number of nodes\t"<<Num_of_Nodes<<endl;
    
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
