#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>




using namespace std::complex_literals;

using namespace std;



vector<double> Membrane::get_ulmYlm_vectorlist_for_mesh(){
    
    //    double radius = sqrt( surface_area_voronoi/(M_PI*4) );
    
    //        calculate_volume_and_surface_area();
    //        double radius = cbrt( 3*volume/(M_PI*4) );
    //
    double radius=0;
    //    double radius = sqrt( surface_area_voronoi/(M_PI*4) );
    for(int i=0;i<Num_of_Nodes;i++){
        radius+=spherical_positions[i][0];
    }
    radius/=Num_of_Nodes;
//    cout<<radius<<endl;
    
    vector<double> radius_vectorlist;
    radius_vectorlist.resize(Num_of_Nodes);
    for(int i=0;i<Num_of_Nodes;i++){
        radius_vectorlist[i] =(spherical_positions[i][0]/radius-1);
    }
    return radius_vectorlist;
}


vector<std::complex<double> > Membrane::get_ylm_vectorlist_for_mesh(int ell, int m, bool complex_conjugate){
    
    vector<std::complex<double> > ylm_vectorlist;
    ylm_vectorlist.resize(Num_of_Nodes);
    for(int i=0;i<Num_of_Nodes;i++){
        ylm_vectorlist[i] = calc_complex_ylmthetaphi(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
        if (complex_conjugate) {
            double im = imag(ylm_vectorlist[i]);
            ylm_vectorlist[i].imag(-im);
        }
    }
    return ylm_vectorlist;
}

vector<double> Membrane::get_real_ylm_vectorlist_for_mesh(int ell, int m){
    
    vector<double> ylm_vectorlist;
    ylm_vectorlist.resize(Num_of_Nodes);
    for(int i=0;i<Num_of_Nodes;i++){
        ylm_vectorlist[i] = calc_real_ylmthetaphi(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
    }
    return ylm_vectorlist;
}



std::complex<double> Membrane::calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<std::complex<double> > vectorlist2){
    std::complex<double> sum=0;
    for(int i=0;i<Num_of_Nodes;i++){
        sum += vectorlist1[i]*vectorlist2[i]*node_voronoi_area[i]
        /(spherical_positions[i][0]*spherical_positions[i][0]);;
    }
    return sum;
}
std::complex<double> Membrane::calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<double> vectorlist2){
    std::complex<double> sum=0;
    for(int i=0;i<Num_of_Nodes;i++){
        sum += vectorlist1[i]*vectorlist2[i]*node_voronoi_area[i]
        /(spherical_positions[i][0]*spherical_positions[i][0]);;
    }
    return sum;
}

double Membrane::calc_vectorlist_vectorlist_surface_integral(vector<double> vectorlist1, vector<double> vectorlist2){
    double sum=0;
    for(int i=0;i<Num_of_Nodes;i++){
        sum += vectorlist1[i]*vectorlist2[i]*node_voronoi_area[i]
        /(spherical_positions[i][0]*spherical_positions[i][0]);;
    }
    return sum;
}
