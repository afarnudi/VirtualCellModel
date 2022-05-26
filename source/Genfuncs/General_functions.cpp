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
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

vector<double> crossvector(vector<double> a, vector<double> b) // cross porduct
{
    if (a.size()!=b.size()) {
        string errorMessage = "crossvector: Both vectors must have dimention 3. I got "+to_string(a.size())+" and "+to_string(b.size())+".";
        throw std::runtime_error(errorMessage);
    }
    vector<double> c(3,0);
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
    return c;
}

vector<double> vector_subtract(vector<double> vec1, vector<double> vec2){
    if (vec1.size()!=vec2.size()) {
        string errorMessage = "vector_subtract: Both vectors must have the same dimention. I got "+to_string(vec1.size())+" and "+to_string(vec2.size())+".";
        throw std::runtime_error(errorMessage);
    }
    vector<double> subtract(vec1.size(),0);
    for (int i=0; i<subtract.size(); i++) {
        subtract[i]=vec1[i]-vec2[i];
    }
    return subtract;
}


double innerproduct(double n1[3],double n2[3])
{
    return n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
}

double innerproduct(vector<double> n1,vector<double> n2)
{
    if (n1.size()!=n2.size()) {
        string errorMessage = "innerproduct: Both vectors must have the same dimention. I got "+to_string(n1.size())+" and "+to_string(n2.size())+".";
        throw std::runtime_error(errorMessage);
    }
    double sum=0;
    for (int i=0; i<n1.size(); i++) {
        sum+=n1[i]*n2[i];
    }
    return sum;
}

double vector_length(double v[3]) // calculate length of vector
{
    return  sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );
}
double vector_length(vector<double> v) // calculate length of vector
{
    double element_sum=0;
    for (int i=0; i<v.size(); i++) {
        element_sum+=v[i]*v[i];
    }
    return  sqrt( element_sum );
}


double vector_length_squared(double v[3]) // calculate length of vector
{
    return  v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}
double vector_length_squared(vector<double> v) // calculate length of vector
{
    double element_sum=0;
    for (int i=0; i<v.size(); i++) {
        element_sum+=v[i]*v[i];
    }
    return  element_sum;
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

vector<double> normalise_vector(vector<double> vec){
    vector<double> normalised_vec(vec.size(),0);
    double vec_length = vector_length(vec);
    for (int i=0; i<vec.size(); i++) {
        normalised_vec[i]=vec[i]/vec_length;
    }
    return normalised_vec;
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
//        cout<<line<<endl;
//        for (auto &i:splitline){
//            cout<<i<<endl;
//        }
        checkpointPlatformName=splitline[3];
        
//        exit(0);
    }
    
    
    
    return reportfilename;
}

double get_double_value(std::string value, std::string parser_name, std::string value_name, std::string example_inputs){
    double test;
    try {
        test = stod(value);
    } catch (...) {
        string errorMessage = TWARN;
        errorMessage+=parser_name + ": Invalid input for the \""+value_name+"\" (";
        errorMessage+=TFILE;
        errorMessage+=value;
        errorMessage+=TWARN;
        errorMessage+="). Please try again.\nExample inputs: "+example_inputs;
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    return test;
}

double get_int_value(std::string value, std::string parser_name, std::string value_name, std::string example_inputs){
    double test;
    try {
        test = stoi(value);
    } catch (...) {
        string errorMessage = TWARN;
        errorMessage+=parser_name + ": Invalid input for the \""+value_name+"\" (";
        errorMessage+=TFILE;
        errorMessage+=value;
        errorMessage+=TWARN;
        errorMessage+="). Integer input required. Please try again.\nExample inputs: "+example_inputs;
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    return test;
}

bool get_bool_value(std::string value, std::string parser_name, std::string value_name){
    bool test;
    if(value=="true"){
        test=true;
    } else if (value=="false"){
        test=false;
    } else {
        string errorMessage = TWARN;
        errorMessage+=parser_name + ": Invalid input for the \""+value_name+"\" (";
        errorMessage+=TFILE;
        errorMessage+=value;
        errorMessage+=TWARN;
        errorMessage+="). Boolian input required. Please try again with \"true\" or \"false\".";
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
    return test;
}


void split(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
    /* Store the original string in the array, so we can loop the rest
     * of the algorithm. */
    tokens.clear();
    tokens.push_back(str);

    // Store the split index in a 'size_t' (unsigned integer) type.
    size_t splitAt;
    // Store the size of what we're splicing out.
    size_t splitLen = splitBy.size();
    // Create a string for temporarily storing the fragment we're processing.
    std::string frag;
    // Loop infinitely - break is internal.
    while(true)
    {
        /* Store the last string in the vector, which is the only logical
         * candidate for processing. */
        frag = tokens.back();
        /* The index where the split is. */
        splitAt = frag.find(splitBy);
        // If we didn't find a new split point...
        if(splitAt == std::string::npos)
        {
            // Break the loop and (implicitly) return.
            break;
        }
        /* Put everything from the left side of the split where the string
         * being processed used to be. */
        tokens.back() = frag.substr(0, splitAt);
        /* Push everything from the right side of the split to the next empty
         * index in the vector. */
        tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
    }
}
