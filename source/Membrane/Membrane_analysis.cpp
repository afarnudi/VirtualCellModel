#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>


using namespace std::complex_literals;

using namespace std;



void Membrane::load_pdb_frame(int frame, int analysis_averaging_option,int znode, int ynode){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] = pdb_frames[frame][i][j];
        }
    }
    z_node_index = znode;
    y_node_index = ynode;
    update_COM_position();
    set_com_to_zero();
    if (analysis_averaging_option == 1) {
        srand (time(NULL));
        double scratch = rand();
        scratch = rand();
    }else if(analysis_averaging_option == 2){
        if(z_node_index == -1 || y_node_index == -1){
            cout<<"No nodes specified for aligning.\n the nodes must be specfied with the '-align_axes flag z y' where z is the node index to be aligned with the  z axis (when the com is in the origin (0,0,0) )and y is the index of the node that will then be rotated to lie on the zy plane.\n";
            exit(EXIT_FAILURE);
        }
        rotate_particle_to_axes();
    }
    
    
    update_spherical_positions();
    calculate_surface_area_with_voronoi();
}



void Membrane::calculate_ulm(int ell_max, int analysis_averaging_option){
    if(ulm_avg.size() != ell_max+1){
        ulm_avg.clear();
        ulm_std.clear();
        ulm_avg.resize(ell_max+1);
        ulm_std.resize(ell_max+1);
        for (int ell=0; ell<ell_max+1; ell++) {
            ulm_avg[ell].resize(2*ell+1,0);
            ulm_std[ell].resize(2*ell+1,0);
        }
        cout<<"cleared ulm\n";
    }
    
    
    if(analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    
    
    double radius = sqrt( surface_area_voronoi/(M_PI*4) );
    //    cout<<"radius voronoi = "<<radius<<endl;
//    double r=0;
//    for(int i=0;i<Num_of_Nodes;i++){
//        r+=spherical_positions[i][0];
//    }
//    r/=Num_of_Nodes;
    //    cout<<"radius spheric = "<<r<<endl;
    
    std::complex<double> ylm;
    std::complex<double> ylm_cc;
    
    double voronoi_to_omega_multiplyer = 1./(radius);
    double f_theta_phi;
    
    
    vector<vector< std::complex<double> > > ulm_avg_frame;
    ulm_avg_frame.resize(ell_max+1);
    for (int ell=0; ell<ell_max+1; ell++) {
        ulm_avg_frame[ell].resize(2*ell+1);
        
        for (int m=0; m<2*ell+1; m++) {
            ulm_avg_frame[ell][m].real(0);
            ulm_avg_frame[ell][m].imag(0);
        }
    }
    
    for (int ell=0; ell<ell_max+1; ell++) {
        //        cout<<ell<<" out of "<<ell_max<<"\r";
        for (int m=-ell; m<ell+1; m++) {
            for(int i=0;i<Num_of_Nodes;i++){
                ylm = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
                ylm_cc = ylm;
                ylm_cc.imag(-1*imag(ylm));

                f_theta_phi = (spherical_positions[i][0]-radius);
                double multi = f_theta_phi*voronoi_to_omega_multiplyer
                                          *node_voronoi_area[i]
                                          /(spherical_positions[i][0]*spherical_positions[i][0]);
                ylm_cc.real(real(ylm) * multi);
                ylm_cc.imag(imag(ylm) * multi);

                ulm_avg_frame[ell][m+ell] += ylm_cc;
            }

        }
    }
    
    
    for (int ell=0; ell<ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            double ulm = abs(ulm_avg_frame[ell][m+ell]);
            ulm_avg[ell][m+ell] += ulm*ulm;
            ulm_std[ell][m+ell] += ulm*ulm*ulm*ulm;
            
        }
    }
    
    //    cout<<endl;
}

void Membrane::write_ulm(int ell_max, string traj_name, double num_frames, string extension){
    
    traj_name.erase(traj_name.end()-4,traj_name.end() );
    string traj_file_name = traj_name + extension;
    
    std::ofstream wdata;
    wdata.open(traj_file_name.c_str(), std::ios::app);
    
    
    
    for (int ell=0; ell<ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            
            wdata<<ulm_avg[ell][m+ell]<<"\t";
        }
        wdata<<"\n";
        
    }
    for (int ell=0; ell<ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            
            wdata<<ulm_std[ell][m+ell]<<"\t";
        }
        wdata<<"\n";
    }
    cout<<traj_name<<" written"<<endl;
}
