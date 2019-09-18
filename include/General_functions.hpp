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
#include "General_constants.h"


void crossvector( double c[3],double d[3],double b[3] ); // cross porduct
double innerproduct(double n1[3],double n2[3]);
double vector_length(double v[3]); // calculate length of vector
double vector_length_squared(double v[3]);
double sign_function(double x);
double periodiccondition(double dx );
void Vector_transformation (double MV[3],double  M[3][3] ,double V[3]);
void matrix_inverse (double mat[3][3]);

//Calculate the A, B, C coefficeints for the surface equation: Ax + By + Cz = -1
void calc_surface_coefficeints (double points[3][3], double &A, double &B, double &C);



void print_wall_clock_time(double sim_duration_per_sec);
void print_real_time(std::chrono::time_point<std::chrono::steady_clock> chrono_clock_start,
                     std::chrono::time_point<std::chrono::steady_clock> chrono_clock_end);
void print_system_time(std::chrono::time_point<std::chrono::system_clock> chrono_clock_start,
                       std::chrono::time_point<std::chrono::system_clock> chrono_clock_end);
#endif /* General_functions_hpp */
