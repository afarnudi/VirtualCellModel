#include <sstream>

#include "Membrane.h"
#include "Configfile.hpp"

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
//    cout<<args.Mesh_files[file_index]<<endl;
    if (args.analysis!="E") {
        analysis_init(args.Mesh_files[file_index]);
    }
    
    
    std::ifstream read_pdb;
    read_pdb.open(args.analysis_filename.c_str());
    if (!read_pdb) {
        cout << "Unable to read "<<args.analysis_filename<<endl;
        exit(EXIT_FAILURE);
    }
    string line;
    read_pdb.seekg(std::ios::beg);
    int num_of_frames = args.framelimits_end-args.framelimits_beg;
    
    analysis_coord_frames_time.clear();
    analysis_coord_frames_time.resize(num_of_frames);
    
    analysis_coord_frames.clear();
    analysis_coord_frames.resize(num_of_frames);
    for (int i=0; i<num_of_frames; i++) {
        analysis_coord_frames[i].resize(Num_of_Nodes);
        for (int j=0; j<Num_of_Nodes; j++) {
            analysis_coord_frames[i][j].resize(3);
            analysis_coord_frames[i][j][0]=0;
            analysis_coord_frames[i][j][1]=0;
            analysis_coord_frames[i][j][2]=0;
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
        
        auto split = split_and_check_for_comments(line, "Import PDB: ");
        string time = split[2];
        time.erase(time.begin(),time.begin()+5);
        analysis_coord_frames_time[frame_index]=stod(time);
        
        int local_atom_count=0;
        for (int node_index= 0; node_index<args.num_atoms_per_frame; node_index++) {
            
            getline(read_pdb, line);
            std::istringstream iss(line);
            vector<string> split(std::istream_iterator<string>{iss}, std::istream_iterator<string>());
            if (split[2]==label) {
                if(split.size()==11){
                    analysis_coord_frames[frame_index][local_atom_count][0] = stod(split[6]);
                    analysis_coord_frames[frame_index][local_atom_count][1] = stod(split[7]);
                    analysis_coord_frames[frame_index][local_atom_count][2] = stod(split[8]);
                } else if(split.size()==10){
                    if (split[6].size()<=9){
                        
                        analysis_coord_frames[frame_index][local_atom_count][0] = stod(split[6]);
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
                        analysis_coord_frames[frame_index][local_atom_count][1] = stod(coor);
                        
                        coor = split[7];
                        coor.erase(coor.begin() + 0,coor.begin()+ it );
                        analysis_coord_frames[frame_index][local_atom_count][2] = stod(coor);
                        
                        
                        
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
                        analysis_coord_frames[frame_index][local_atom_count][0] = stod(coor);
                        
                        coor = split[6];
                        coor.erase(coor.begin() + 0,coor.begin()+ it );
                        analysis_coord_frames[frame_index][local_atom_count][1] = stod(coor);
                        analysis_coord_frames[frame_index][local_atom_count][2] = stod(split[7]);
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
                    analysis_coord_frames[frame_index][local_atom_count][0] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + it[1], coor.end());
                    coor.erase(coor.begin() + 0,coor.begin()+ it[0] );
                    analysis_coord_frames[frame_index][local_atom_count][1] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + 0,coor.begin()+ it[1] );
                    analysis_coord_frames[frame_index][local_atom_count][2] = stod(coor);
                    
                }
                local_atom_count++;
            }
            

        }
        if (local_atom_count != Num_of_Nodes) {
            cout<<"Number of atoms in frame, "<<frame_index+args.framelimits_beg<<", does not match with the number of atoms in previous frames. Expected "<<Num_of_Nodes<<" got "<<local_atom_count<<endl;
            exit(EXIT_FAILURE);
        }
        getline(read_pdb, line);
    }
}

void Membrane::import_xyz_frames(ArgStruct_Analysis args, int file_index){
    string label = args.membane_labels[file_index];
    Num_of_Nodes      = get_xyz_num_of_atoms(args.analysis_filename, label);
    
    if (label == "") {
        label = get_xyz_first_label(args.analysis_filename);
    }
    if (!generalParameters.Testmode) {
        cout<<Num_of_Nodes<<" "<<label<<" nodes in "<<args.analysis_filename<<endl;
    }
//    cout<<args.Mesh_files[file_index]<<endl;
    if (args.analysis!="E") {
        analysis_init(args.Mesh_files[file_index]);
    }
    
    
    std::ifstream read_xyz;
    read_xyz.open(args.analysis_filename.c_str());
    if (!read_xyz) {
        cout << "Unable to read "<<args.analysis_filename<<endl;
        exit(EXIT_FAILURE);
    }
    string line;
    read_xyz.seekg(std::ios::beg);
    int num_of_frames = args.framelimits_end-args.framelimits_beg;
    
    analysis_coord_frames_time.clear();
    analysis_coord_frames_time.resize(num_of_frames);
    
    analysis_coord_frames.clear();
    analysis_coord_frames.resize(num_of_frames);
    for (int i=0; i<num_of_frames; i++) {
        analysis_coord_frames[i].resize(Num_of_Nodes);
        for (int j=0; j<Num_of_Nodes; j++) {
            analysis_coord_frames[i][j].resize(3);
            analysis_coord_frames[i][j][0]=0;
            analysis_coord_frames[i][j][1]=0;
            analysis_coord_frames[i][j][2]=0;
        }
    }
    //Skip frames that are not needed
    for (int frame_index= 0; frame_index<args.framelimits_beg; frame_index++){
        getline(read_xyz, line);
        getline(read_xyz, line);
        for (int node_index= 0; node_index<args.num_atoms_per_frame; node_index++) {
            getline(read_xyz, line);
        }
    }
    
    for (int frame_index= 0; frame_index<num_of_frames; frame_index++){

        getline(read_xyz, line);
        getline(read_xyz, line);
        
        auto split = split_and_check_for_comments(line, "Import XYZ: ");
        string time = split[1];
        analysis_coord_frames_time[frame_index]=stod(time);

        int local_atom_count=0;
        double readCoord[3];
        string readLable;
        for (int node_index= 0; node_index<args.num_atoms_per_frame; node_index++) {
            read_xyz>>readLable;
            read_xyz>>readCoord[0]>>readCoord[1]>>readCoord[2];
            if (readLable==label) {
                analysis_coord_frames[frame_index][local_atom_count][0]=readCoord[0];
                analysis_coord_frames[frame_index][local_atom_count][1]=readCoord[1];
                analysis_coord_frames[frame_index][local_atom_count][2]=readCoord[2];
                local_atom_count++;
            }
        }
        if (local_atom_count != Num_of_Nodes) {
            cout<<"Number of atoms in frame, "<<frame_index+args.framelimits_beg<<", does not match with the number of atoms in previous frames. Expected "<<Num_of_Nodes<<" got "<<local_atom_count<<endl;
            exit(EXIT_FAILURE);
        }
    }
}
