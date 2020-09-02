#include "General_functions.hpp"
#include <iostream>
#include <sstream>
#include <string>

int get_pdb_num_of_atoms(std::string filename, std::string label){
    std::ifstream read_pdb;
    read_pdb.open(filename.c_str());
    read_pdb.seekg(std::ios::beg);
    if (!read_pdb.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    int num_atoms = 0;
    bool finished_frame=false;
    getline(read_pdb, line);
    getline(read_pdb, line);
    while(!finished_frame){
        getline(read_pdb, line);
        std::istringstream iss(line);
        std::vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
        if (split[2]==label) {
            num_atoms++;
        } else {
            finished_frame=true;
        }
        
    }
    read_pdb.close();
    
    return num_atoms;
}

int get_pdb_num_of_atoms(std::string filename){
    std::ifstream read_pdb;
    read_pdb.open(filename.c_str());
    read_pdb.seekg(std::ios::beg);
    if (!read_pdb.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    
    int num_atoms = 0;
    bool finished_frame=false;
    getline(read_pdb, line);
    getline(read_pdb, line);
    while(!finished_frame){
        getline(read_pdb, line);
        std::istringstream iss(line);
        std::vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
        if (split[0]=="ATOM") {
            num_atoms++;
        } else {
            finished_frame=true;
        }
        
    }
    read_pdb.close();
    
    return num_atoms;
}

int get_pdb_num_of_frames(std::string filename, int num_atoms){
    std::ifstream read_pdb;
    read_pdb.open(filename.c_str());
    read_pdb.seekg(std::ios::beg);
    if (!read_pdb.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    
    int num_frames = 0;
    while(getline(read_pdb, line)){
        num_frames++;
    }
    read_pdb.close();
//    std::cout<<"num_atoms = "<<num_atoms<<std::endl;
//    std::cout<<"num_lines = "<<num_frames<<std::endl;
    return num_frames/(num_atoms+3);
}



std::string get_pdb_first_label(std::string filename){
    std::ifstream read_pdb;
    read_pdb.open(filename.c_str());
    read_pdb.seekg(std::ios::beg);
    
    std::string line;
    std::string label;
    
    getline(read_pdb, line);
    getline(read_pdb, line);
    getline(read_pdb, line);
    
    std::istringstream iss(line);
    std::vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    label= split[2];
    read_pdb.close();
    
    return label;
}


