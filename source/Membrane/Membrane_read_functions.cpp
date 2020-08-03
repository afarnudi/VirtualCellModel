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
    shift_position(Shift_in_X_direction, Shift_in_Y_direction, Shift_in_Z_direction);
    cout<<"\n\nMembrane class initiated.\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}


void Membrane::read_ply_file (std::string ply_file)
{
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(ply_file.c_str()); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary ply file intigers. We never use these intigers in the actual programme.
    std::string temp_string;
    for (int i=0; i<18; i++)
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
    for (int i=0; i<11; i++)
    {
        read>> temp_string;
    }
    read>>Num_of_Triangles;
    
    for (int i=0; i<6; i++)
    {
        read>> temp_string;
    }
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
    if(fixing_com){
        read>>temp_string;
        read>>Num_of_Free_Bonds;
        vector<int> extra;
        extra.resize(2);
        for(int i=0; i<Num_of_Free_Bonds; i++){
        read>>extra[0];
        read>>extra[1];
        Free_Bonds.push_back(extra);
        }
    cout<<"Num_of_Free_Bonds  "<<Num_of_Free_Bonds<<endl;
    }
    
}


#include <sstream>
using std::string;

int Membrane::import_pdb_frames(string filename){
    
    
    int num_of_frames = count_pdb_frames(filename, Num_of_Nodes);
//    cout<<"I'm out\n";
    ifstream read_pdb;
    read_pdb.open(filename.c_str());
    if (!read_pdb) {
        cout << "Unable to read "<<filename<<endl;
        exit(EXIT_FAILURE);
    }
    string line;
    read_pdb.seekg(std::ios::beg);
    
    if (!GenConst::Testmode) {
        cout<<"num of frames = "<<num_of_frames<<endl;
    }
    
    pdb_frames.resize(num_of_frames);
    for (int i=0; i<num_of_frames; i++) {
        pdb_frames[i].resize(Num_of_Nodes);
        for (int j=0; j<Num_of_Nodes; j++) {
            pdb_frames[i][j].resize(3);
            pdb_frames[i][j][0]=0;
            pdb_frames[i][j][1]=0;
            pdb_frames[i][j][2]=0;
        }
    }

    for (int frame_index= 0; frame_index<num_of_frames; frame_index++){

        getline(read_pdb, line);
        getline(read_pdb, line);

        for (int node_index= 0; node_index<Num_of_Nodes; node_index++) {
            
            getline(read_pdb, line);
            std::istringstream iss(line);
            vector<string> split(std::istream_iterator<string>{iss}, std::istream_iterator<string>());
            
            if(split.size()==11){
                pdb_frames[frame_index][node_index][0] = stod(split[6]);
                pdb_frames[frame_index][node_index][1] = stod(split[7]);
                pdb_frames[frame_index][node_index][2] = stod(split[8]);
            } else if(split.size()==10){
                if (split[6].size()<=9){
                    
                    pdb_frames[frame_index][node_index][0] = stod(split[6]);
                    int it=-1;
                    for(int i=0; i<split[7].size();i++){
                        if (split[7][i] == '.' ){
                            it=i+4;
                            break;
                        }
                    }
                    
                    
                    string coor;
                    coor = split[7];
                    coor.erase(coor.begin() + it, coor.end());
                    pdb_frames[frame_index][node_index][1] = stod(coor);
                    
                    coor = split[7];
                    coor.erase(coor.begin() + 0,coor.begin()+ it );
                    pdb_frames[frame_index][node_index][2] = stod(coor);
                    
                    
                    
                } else {
                    int it=-1;
                    for(int i=0; i<split[6].size();i++){
                        if (split[6][i] == '.' ){
                            it=i+4;
                            break;
                        }
                    }
                    
                    string coor;
                    coor = split[6];
                    coor.erase(coor.begin() + it, coor.end());
                    pdb_frames[frame_index][node_index][0] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + 0,coor.begin()+ it );
                    pdb_frames[frame_index][node_index][1] = stod(coor);
                    pdb_frames[frame_index][node_index][2] = stod(split[7]);
                }
            } else if(split.size()==9){
                int it[3] = {-1,-1,-1};
                int it_index=0;
                for(int i=0; i<split[6].size();i++){
                    if (split[6][i] == '.' ){
                        it[it_index]=i+4;
                        it_index++;
                    }
                }
                
                string coor;
                coor = split[6];
                coor.erase(coor.begin() + it[0], coor.end());
                pdb_frames[frame_index][node_index][0] = stod(coor);
                
                coor = split[6];
                coor.erase(coor.begin() + it[1], coor.end());
                coor.erase(coor.begin() + 0,coor.begin()+ it[0] );
                pdb_frames[frame_index][node_index][1] = stod(coor);
                
                coor = split[6];
                coor.erase(coor.begin() + 0,coor.begin()+ it[1] );
                pdb_frames[frame_index][node_index][2] = stod(coor);
                
            }

        }
        getline(read_pdb, line);
    }
    return num_of_frames;
}
