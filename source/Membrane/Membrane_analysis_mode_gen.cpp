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

#include <boost/math/special_functions/legendre.hpp>
#include <cmath>
double Membrane::calc_assoc_legendre(int ell, int m, double x){
    return boost::math::legendre_p(ell, m, x);
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

    }
    set_com_to_zero();
    calculate_surface_area_with_voronoi();
}

void Membrane::add_ulm_mode_real(int Ell, int M, double uLM, double radius){
    
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
        vector<double> add_coordinates;
        add_coordinates.resize(3,0);
        spherical_positions[i][0] += radius*(uLM*ylmlist[i]);
        
        Node_Position[i] = convert_spherical_to_cartesian(spherical_positions[i][0], spherical_positions[i][1], spherical_positions[i][2]);
        
    }
    set_com_to_zero();
    calculate_surface_area_with_voronoi();
}




void Membrane::randomundulationgenerator(void){
    
    
    
    
    for (int i=0; i<NumberOfRandomModes; i++) {
        double r = ((double) rand() / (RAND_MAX));
        int Ell = rand() % (EllMaxOfRandomModes+1);
        int M   = rand() % (Ell+1);
        if (r<0.5) {
            M*=-1;
        }
        update_spherical_positions();
        calculate_dOmega();
        calculate_surface_area_with_voronoi();
        add_ulm_mode_real(Ell, M, UlmOfRandomModes, Radius);
        cout<<"Ulm "<<UlmOfRandomModes<<"  Ell "<<Ell<<"  M "<<M<<endl;
    }
    
    
}

void Membrane::generate_undulations(void){
    update_spherical_positions();
    calculate_dOmega();
    calculate_surface_area_with_voronoi();
    add_ulm_mode_real(EllOfGeneratedMode, MOfGeneratedMode, AmplitudeOfGeneratedMode, Radius);
}

void Membrane::randomundulationgenerator(int numberOfRandomModes,
                                         int ellMaxOfRandomModes,
                                         double ulmOfRandomModes){
    for (int i=0; i<numberOfRandomModes; i++) {
        double r = ((double) rand() / (RAND_MAX));
        int Ell = 2 + rand() % (ellMaxOfRandomModes - 1);
        int M   = rand() % (Ell+1);
        if (rand()<0.5) {
            M*=-1;
        }
        update_spherical_positions();
        calculate_dOmega();
        calculate_surface_area_with_voronoi();
        add_ulm_mode_real(Ell, M, ulmOfRandomModes, Radius);
        cout<<"Ulm "<<ulmOfRandomModes<<"  Ell "<<Ell<<"  M "<<M<<endl;
    }
    
    
}
