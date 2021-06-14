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
    double p1p2[3], p1p3[3], n[3];//, D;
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

std::vector<double> convert_cartesian_to_spherical(std::vector<double> xyz){
    std::vector<double> r_theta_phi;
    r_theta_phi.resize(3,0);
    
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > p1;
    bg::model::point<double, 3, bg::cs::cartesian> p3(xyz[0],xyz[1],xyz[2]);
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

std::vector<double> convert_spherical_to_cartesian(std::vector<double> r_theta_phi){
    std::vector<double> xyz;
    xyz.resize(3,0);
    
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > p1(r_theta_phi[2],r_theta_phi[1],r_theta_phi[0]);
    bg::model::point<double, 3, bg::cs::cartesian> p3;
    bg::transform(p1, p3);
    
    xyz[0] = p3.get<0>();
    xyz[1] = p3.get<1>();
    xyz[2] = p3.get<2>();
    return xyz;
}


string check_if_file_exists(string filename){
    ifstream infile(filename.c_str());
    bool loop=infile.is_open();
    string entree = filename;
    while( !loop ){
        cout<<"'"<<TBOLD<<TWARN<<filename<<TRESET<<"' does not exist. Please try again:\n"<<TBLINK<<">> "<<TRESET;
        cin>>entree;
        ifstream infiletemp(entree.c_str());
        loop=infiletemp.is_open();
    }
    return entree;
}


#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

std::vector<std::string> splitstring(std::string line, char delimiter){
    std::stringstream ss(line);
    std::string s;
    const char delim = delimiter;
    std::vector<std::string> splitpath;
    while (std::getline(ss, s, delim)) {
        splitpath.push_back(s);
    }
    return splitpath;
}

#include <boost/filesystem.hpp>
std::string find_resume_config(std::string resumePath, std::string &checkpointPath, std::string &checkpointPlatformName){
    std::vector<std::string> splitpath = splitstring(resumePath, '/');
    if(splitpath.size()==0){
        string errorMessage = TWARN;
        errorMessage+="error parsing options: argument -r: Please provide a directory path.";
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    
    string reportfilename, hardwarereport;
    for (int i=0; i<splitpath.size(); i++) {
        reportfilename+=splitpath[i]+"/";
        checkpointPath+=splitpath[i]+"/";
        hardwarereport+=splitpath[i]+"/";
        generalParameters.trajectory_file_name+=splitpath[i]+"/";
        generalParameters.buffer_file_name+=splitpath[i]+"/";
    }
    reportfilename+=splitpath[splitpath.size()-1]+"_report.txt";
    checkpointPath+=splitpath[splitpath.size()-1]+"_Checkpoint";
    hardwarereport+=splitpath[splitpath.size()-1]+"_hardware_runtime.txt";
    generalParameters.trajectory_file_name+=splitpath[splitpath.size()-1];
    generalParameters.buffer_file_name+="buffs/"+splitpath[splitpath.size()-1];
    
    ifstream reportfile(reportfilename.c_str());
    bool fileExists=reportfile.is_open();
    if( !fileExists ){
        string errorMessage = TWARN;
        errorMessage+="error parsing options: resume: report file not found:\n";
        errorMessage+=TFILE;
        errorMessage+=reportfilename;
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    ifstream checkpointfile(checkpointPath.c_str());
    fileExists=checkpointfile.is_open();
    if( !fileExists ){
        string errorMessage = TWARN;
        errorMessage+="error parsing options: resume: checkpoint file not found:\n";
        errorMessage+=TFILE;
        errorMessage+=checkpointPath;
        errorMessage+= TRESET;
        cout<< errorMessage<<endl;
        
        cout<<"Will try backup.."<<endl;
        string backupCheckpointPath = checkpointPath+"Backup";
        ifstream backupcheckpointfile(backupCheckpointPath.c_str());
        fileExists=backupcheckpointfile.is_open();
        if( !fileExists ){
            string errorMessage = TWARN;
            errorMessage+="Both the checkpoint and the backup are missing (in the simulation directory). This simulation cannot be resumed.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
    }
    
    ifstream hardwarefile(hardwarereport.c_str());
    fileExists=hardwarefile.is_open();
    if( !fileExists ){
        cout<<TWARN<<"warning hardware report file not in the simulation path. Make sure the simulation is resuming on the same machine.\n";
        cout<<hardwarereport<<endl<<TRESET;
    } else {
        string line;
        getline (hardwarefile,line);
        string machineName = exec("hostname");
        
        vector<string> splitline = splitstring(machineName, '\n');
        machineName=splitline[0];
        if (machineName!=line) {
            
            string hardwarereportbackup = hardwarereport+"Backup";
            ifstream hardwarefilebackup(hardwarereport.c_str());
            getline (hardwarefilebackup,line);
            if (machineName!=line) {
                string errorMessage = TWARN;
                errorMessage+="error parsing options: resume: Resume machine name does not match:\n";
                errorMessage+="Simulation was run on :    ";
                errorMessage+=TFILE;
                errorMessage+=line+"\n";
                errorMessage+=TWARN;
                errorMessage+="Simulation is resuming on: ";
                errorMessage+=TFILE;
                errorMessage+=machineName+"\n";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            } else {
                boost::filesystem::copy(hardwarereportbackup, hardwarereport);
            }
            
            
        }
        getline (hardwarefile,line);
        splitline = splitstring(line, ' ');
        cout<<line<<endl;
        for (auto &i:splitline){
            cout<<i<<endl;
        }
        checkpointPlatformName=splitline[3];
        
//        exit(0);
    }
    
    
    
    return reportfilename;
}
