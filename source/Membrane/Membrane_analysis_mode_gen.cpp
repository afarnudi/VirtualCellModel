#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
using namespace std::complex_literals;
using namespace std;



std::complex<double> Membrane::calc_complex_ylmthetaphi(int l,  int  m, double theta, double phi){
    return boost::math::spherical_harmonic(l,m,theta,phi);
}


double Membrane::calc_real_ylmthetaphi(int l,  int  m, double theta, double phi){
    std::complex<double> ylm;
    
    double realYLM=0;
    
    if (m<0){
        ylm = calc_complex_ylmthetaphi(l,-m,theta,phi);
        
        realYLM = imag(ylm)*pow(-1,m)*sqrt(2);
        
    } else if (m==0){
        ylm = calc_complex_ylmthetaphi(l,0,theta,phi);
        
        realYLM = real(ylm);
    } else {
        
        ylm = calc_complex_ylmthetaphi(l,m,theta,phi);
        
        realYLM = real(ylm)*pow(-1,m)*sqrt(2);
        
    }
    return realYLM;
}




void Membrane::generate_ulm_mode(int Ell, int M, double uLM, double radius){
    if (Ell<0){
        cout<<"\nl must be positive"<<endl;
        exit(EXIT_FAILURE);
    }
    if (abs(M)>Ell){
        cout<<"\n |m| cannot be larger than l"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<std::complex<double> > ylmlist = get_ylm_vectorlist_for_mesh(Ell, M, false);
    
    for(int i=0;i<Num_of_Nodes;i++){
        spherical_positions[i][0] = radius*(1+uLM*real(ylmlist[i]));
        Node_Position[i]=convert_spherical_to_cartesian(spherical_positions[i][0], spherical_positions[i][1], spherical_positions[i][2]);
//        Node_Position[i][0]=spherical_positions[i][0]*sin(spherical_positions[i][1])
//        *cos(spherical_positions[i][2]);
//        Node_Position[i][1]=spherical_positions[i][0]*sin(spherical_positions[i][1])
//        *sin(spherical_positions[i][2]);
//        Node_Position[i][2]=spherical_positions[i][0]*cos(spherical_positions[i][1]);
        
    }
    
}


void Membrane::generate_ulm_mode_real(int Ell, int M, double uLM, double radius){
    
    if (Ell<0){
        cout<<"\nl must be positive"<<endl;
        exit(EXIT_FAILURE);
    }
    if (abs(M)>Ell){
        cout<<"\n |m| cannot be larger than l"<<endl;
        exit(EXIT_FAILURE);
    }
    
    
    vector<double> ylmlist = get_real_ylm_vectorlist_for_mesh(Ell, M);
    
    for(int i=0;i<Num_of_Nodes;i++){
        
        spherical_positions[i][0] = radius*(1+uLM*ylmlist[i]);
        Node_Position[i]=convert_spherical_to_cartesian(spherical_positions[i][0], spherical_positions[i][1], spherical_positions[i][2]);
//        Node_Position[i][0]=spherical_positions[i][0]*sin(spherical_positions[i][1])
//        *cos(spherical_positions[i][2]);
//        Node_Position[i][1]=spherical_positions[i][0]*sin(spherical_positions[i][1])
//        *sin(spherical_positions[i][2]);
//        Node_Position[i][2]=spherical_positions[i][0]*cos(spherical_positions[i][1]);
        
    }
    set_com_to_zero();
    calculate_surface_area_with_voronoi();
}



