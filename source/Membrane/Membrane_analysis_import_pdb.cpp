#include <sstream>

#include "Membrane.h"

using std::string;

void Membrane::import_pdb_frames(ArgStruct_Analysis args, int file_index){
    string label = args.membane_labels[file_index];
    Num_of_Nodes      = get_pdb_num_of_atoms(args.analysis_filename, label);
    if (label == "") {
        label = get_pdb_first_label(args.analysis_filename);
    }
    if (!generalParameters.Testmode) {
        cout<<Num_of_Nodes<<" "<<label<<" nodes in "<<args.analysis_filename<<endl;
    }
    
//    if (args.analysis_dim==3) {
        analysis_init(args.Mesh_files[file_index]);
//    }
    
    
    std::ifstream read_pdb;
    read_pdb.open(args.analysis_filename.c_str());
    if (!read_pdb) {
        cout << "Unable to read "<<args.analysis_filename<<endl;
        exit(EXIT_FAILURE);
    }
    string line;
    read_pdb.seekg(std::ios::beg);
    int num_of_frames = args.framelimits_end-args.framelimits_beg;
    
    pdb_frames.clear();
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
    //Skip frames that are not needed
    for (int frame_index= 0; frame_index<args.framelimits_beg; frame_index++){
        getline(read_pdb, line);
        getline(read_pdb, line);
        for (int node_index= 0; node_index<args.num_atoms_per_frame; node_index++) {
            getline(read_pdb, line);
        }
        getline(read_pdb, line);
    }
    
    for (int frame_index= 0; frame_index<num_of_frames; frame_index++){

        getline(read_pdb, line);
        getline(read_pdb, line);
        int local_atom_count=0;
        for (int node_index= 0; node_index<args.num_atoms_per_frame; node_index++) {
            
            getline(read_pdb, line);
            std::istringstream iss(line);
            vector<string> split(std::istream_iterator<string>{iss}, std::istream_iterator<string>());
            if (split[2]==label) {
                if(split.size()==11){
                    pdb_frames[frame_index][local_atom_count][0] = stod(split[6]);
                    pdb_frames[frame_index][local_atom_count][1] = stod(split[7]);
                    pdb_frames[frame_index][local_atom_count][2] = stod(split[8]);
                } else if(split.size()==10){
                    if (split[6].size()<=9){
                        
                        pdb_frames[frame_index][local_atom_count][0] = stod(split[6]);
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
                        pdb_frames[frame_index][local_atom_count][1] = stod(coor);
                        
                        coor = split[7];
                        coor.erase(coor.begin() + 0,coor.begin()+ it );
                        pdb_frames[frame_index][local_atom_count][2] = stod(coor);
                        
                        
                        
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
                        pdb_frames[frame_index][local_atom_count][0] = stod(coor);
                        
                        coor = split[6];
                        coor.erase(coor.begin() + 0,coor.begin()+ it );
                        pdb_frames[frame_index][local_atom_count][1] = stod(coor);
                        pdb_frames[frame_index][local_atom_count][2] = stod(split[7]);
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
                    pdb_frames[frame_index][local_atom_count][0] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + it[1], coor.end());
                    coor.erase(coor.begin() + 0,coor.begin()+ it[0] );
                    pdb_frames[frame_index][local_atom_count][1] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + 0,coor.begin()+ it[1] );
                    pdb_frames[frame_index][local_atom_count][2] = stod(coor);
                    
                }
                local_atom_count++;
            }
            

        }
        if (local_atom_count != Num_of_Nodes) {
            cout<<"Number of atoms in the "<<frame_index+args.framelimits_beg<<" does not match with the number of atoms in previous frames. Expected "<<Num_of_Nodes<<" got "<<local_atom_count<<endl;
            exit(EXIT_FAILURE);
        }
        getline(read_pdb, line);
    }
}
