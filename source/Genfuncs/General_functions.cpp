//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "General_functions.hpp"
#include <math.h>
#include <iostream>

void crossvector( double c[3], double a[3],double b[3] ) // cross porduct
{
    c[0]=a[1]*b[2]-a[2]*b[1];    // normal vector to plane (not unitary length)
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

double innerproduct(double n1[3],double n2[3])
{
    return n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
}

double vector_length(double v[3]) // calculate length of vector
{
    return  sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );
}

double vector_length_squared(double v[3]) // calculate length of vector
{
    return  v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

double sign_function(double x){
    if(x<0){
        return -1.0;
    }
    
    return 1.0;
}

void Vector_transformation (double MV[3],double  M[3][3] ,double V[3]){
    MV[0]= M[0][0] * V [0] + M[0][1] * V [1]  +M[0][2] * V [2] ;
    MV[1]= M[1][0] * V [0] + M[1][1] * V [1]  +M[1][2] * V [2] ;
    MV[2]= M[2][0] * V [0] + M[2][1] * V [1]  +M[2][2] * V [2] ;
}


void matrix_inverse (double mat[3][3]){
    int i, j;
    float determinant = 0;
    
    //finding determinant
    for(i = 0; i < 3; i++)
        determinant += (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
    
    double mat_temp[3][3]={0};
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            mat_temp[i][j]=((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant;
        
    }
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            mat[i][j]=mat_temp[i][j];
        
    }
}

void calc_surface_coefficeints_2 (double points[3][3], double &A, double &B, double &C){
    double p1p2[3], p1p3[3], n[3], D;
    for(int i=0; i<3; i++){
        p1p2[i]=points[0][i]-points[1][i];
        p1p3[i]=points[0][i]-points[2][i];
    }
    
    crossvector( n,p1p2,p1p3);
    A=n[0];
    B=n[1];
    C=n[2];
    /*
     D= A*points[0][0]+B*points[0][1]+C*points[0][2];
     if(D!=0){
     A=-A/D;
     B=-B/D;
     C=-C/D;
     }
     check_1= A*points[0][0]+B*points[0][1]+C*points[0][2];
     check_2=A*points[1][0]+B*points[1][1]+C*points[1][2];
     check3=A*points[2][0]+B*points[2][1]+C*points[2][2];
     */
}


#include <boost/geometry.hpp>
namespace bg = boost::geometry;

std::vector<double> convert_cartesian_to_spherical(double x, double y, double z){
    std::vector<double> r_theta_phi;
    r_theta_phi.resize(3,0);
    
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > p1;
    bg::model::point<double, 3, bg::cs::cartesian> p3(x,y,z);
    bg::transform(p3, p1);
    
    r_theta_phi[0] = p1.get<2>();
    r_theta_phi[1] = p1.get<1>();
    r_theta_phi[2] = p1.get<0>();
    return r_theta_phi;
}

std::vector<double> convert_spherical_to_cartesian(double r, double theta, double phi){
    std::vector<double> xyz;
    xyz.resize(3,0);
    
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > p1(phi,theta,r);
    bg::model::point<double, 3, bg::cs::cartesian> p3;
    bg::transform(p1, p3);
    
    xyz[0] = p3.get<0>();
    xyz[1] = p3.get<1>();
    xyz[2] = p3.get<2>();
    return xyz;
}




int count_pdb_frames(std::string filename, int num_atoms){
    std::ifstream read_pdb;
    read_pdb.open(filename.c_str());
    read_pdb.seekg(std::ios::beg);
    if (read_pdb) {
        if (!GenConst::Testmode) {
            std::cout << filename<<" opend successfully. \n\n";
        }
        
    }else{
        std::cout << "Unable to read "<<filename<<std::endl;
    }
    std::string line;
    
    int num_frames = 0;
    while(getline(read_pdb, line)){
        num_frames++;
    }
    read_pdb.close();
    return num_frames/(num_atoms+3);
}

