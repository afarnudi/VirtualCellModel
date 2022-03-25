//
//  General_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
/// \file

#ifndef General_functions_hpp
#define General_functions_hpp
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <iterator>
#include "OpenMM_structs.h"
#include "General_constants.h"


void crossvector( double c[3],double d[3],double b[3] ); // cross porduct
vector<double> crossvector( vector<double> a,vector<double> b); // cross porduct
vector<double> normalise_vector(vector<double> vec);
double innerproduct(double n1[3],double n2[3]);
double innerproduct(vector<double> n1,vector<double> n2);
vector<double> vector_subtract(vector<double> vec1, vector<double> vec2);
double vector_length(double v[3]); // calculate length of vector
double vector_length(vector<double> v);
double vector_length_squared(double v[3]);
double vector_length_squared(vector<double> v);
vector<double> vector_subtract(vector<double> vec1, vector<double> vec2);
double sign_function(double x);
double periodiccondition(double dx );
void Vector_transformation (double MV[3],double  M[3][3] ,double V[3]);
void matrix_inverse (double mat[3][3]);

//Calculate the A, B, C coefficeints for the surface equation: Ax + By + Cz = -1
void calc_surface_coefficeints (double points[3][3], double &A, double &B, double &C);
void calc_surface_coefficeints_2 (double points[3][3], double &A, double &B, double &C);


void print_time(std::string filepath,
                std::string fileBackupPath,
                std::string hardwareReportHeader,
                bool printToScreen,
                clock_t tStart,
                std::chrono::time_point<std::chrono::steady_clock> chrono_steady_clock_start,
                std::chrono::time_point<std::chrono::system_clock> chrono_system_clock_start
                );
void print_wall_clock_time(double sim_duration_per_sec, bool printToScreen);
void print_real_time(std::chrono::time_point<std::chrono::steady_clock> chrono_clock_start,
                     std::chrono::time_point<std::chrono::steady_clock> chrono_clock_end, bool printToScreen);
void print_system_time(std::chrono::time_point<std::chrono::system_clock> chrono_clock_start,
                       std::chrono::time_point<std::chrono::system_clock> chrono_clock_end, bool printToScreen);

std::vector<double> convert_cartesian_to_spherical(double x, double y, double z);
std::vector<double> convert_cartesian_to_spherical(std::vector<double> xyz);
std::vector<double> convert_spherical_to_cartesian(double r, double theta, double phi);
std::vector<double> convert_spherical_to_cartesian(std::vector<double> r_theta_phi);

std::string get_pdb_first_label(std::string filename);
int get_pdb_num_of_atoms(std::string filename, std::string label);
int get_pdb_num_of_atoms(std::string filename);
int get_pdb_num_of_frames(std::string filename, int num_atoms);

std::string get_xyz_first_label(std::string filename);
int get_xyz_num_of_atoms(std::string filename, std::string label);
int get_xyz_num_of_atoms(std::string filename);
int get_xyz_num_of_frames(std::string filename, int num_atoms);

std::string check_if_file_exists(std::string filename);

void assign_project_directories(char* buffer,
                                std::string projectName,
                                std::string &trajPath,
                                std::string &buffPath);

std::string exec(const char* cmd);
std::vector<std::string> splitstring(std::string, char delimiter);

std::string find_resume_config(std::string resumePath, std::string &checkpointPath, std::string &Checkpoint_platformName);

void loadCheckpoint(MyOpenMMData* omm, std::string ckeckpoint_name, std::string ckeckpointBackup_name, bool &usingBackupCheckpoint);
void saveCheckpoint(MyOpenMMData* omm, std::string ckeckpoint_name, std::string ckeckpointBackup_name, bool &usingBackupCheckpoint);

double get_double_value(std::string value,
                        std::string parser_name,
                        std::string value_name,
                        std::string example_inputs);

double get_int_value(std::string value,
                        std::string parser_name,
                        std::string value_name,
                        std::string example_inputs);

bool get_bool_value(std::string value,
                    std::string parser_name,
                    std::string value_name);

void split(std::string str, std::string splitBy, std::vector<std::string>& tokens);

#endif /* General_functions_hpp */
