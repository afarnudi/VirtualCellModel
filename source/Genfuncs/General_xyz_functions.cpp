#include "General_functions.hpp"
#include <iostream>
#include <sstream>
#include <string>

int get_xyz_num_of_atoms(std::string filename, std::string label){
    std::ifstream read_xyz;
    read_xyz.open(filename.c_str());
    read_xyz.seekg(std::ios::beg);
    if (!read_xyz.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    int num_atoms = 0;
    bool finished_frame=false;
    getline(read_xyz, line);
    getline(read_xyz, line);
    while(!finished_frame){
        getline(read_xyz, line);
        if (read_xyz.eof()){
            break;
        }
        std::istringstream iss(line);
        std::vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
        if (split[0]==label) {
            num_atoms++;
        } else {
            finished_frame=true;
        }
        
    }
    read_xyz.close();
    return num_atoms;
}

int get_xyz_num_of_atoms(std::string filename){
    std::ifstream read_xyz;
    read_xyz.open(filename.c_str());
    read_xyz.seekg(std::ios::beg);
    if (!read_xyz.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    getline(read_xyz, line);
    
    int num_atoms = stoi(line);
    read_xyz.close();

    return num_atoms;
}

int get_xyz_num_of_frames(std::string filename, int num_atoms){
    std::ifstream read_xyz;
    read_xyz.open(filename.c_str());
    read_xyz.seekg(std::ios::beg);
    if (!read_xyz.is_open()) {
        std::cout << "Unable to read "<<filename<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;

    int num_frames = 0;
    while(getline(read_xyz, line)){
        num_frames++;
    }
    read_xyz.close();
    return num_frames/(num_atoms+2);
}



std::string get_xyz_first_label(std::string filename){
    std::ifstream read_xyz;
    read_xyz.open(filename.c_str());
    read_xyz.seekg(std::ios::beg);

    std::string line;
    std::string label;

    getline(read_xyz, line);
    getline(read_xyz, line);
    getline(read_xyz, line);

    std::istringstream iss(line);
    std::vector<std::string> split(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    label= split[0];
    read_xyz.close();

    return label;
}


