#include "Membrane.h"
#include "General_functions.hpp"


#include <cstdlib>
#include <complex>




using namespace std::complex_literals;

using namespace std;



void Membrane::load_pdb_frame(int frame, ArgStruct args){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] = pdb_frames[frame][i][j];
        }
    }
    
    update_COM_position();
    set_com_to_zero();
    
    if (args.analysis_averaging_option == 1) {
        srand (time(NULL));
        double scratch = rand();
        scratch = rand();
    }else if(args.analysis_averaging_option == 2){
        if(args.z_node == -1 || args.zy_node == -1){
            cout<<"No nodes specified for aligning.\n the nodes must be specfied with the '-align_axes flag z y' where z is the node index to be aligned with the  z axis (when the com is in the origin (0,0,0) )and y is the index of the node that will then be rotated to lie on the zy plane.\n";
            exit(EXIT_FAILURE);
        }
        rotate_particle_to_axes(args);
    }
    
    
    update_spherical_positions();
    calculate_surface_area_with_voronoi();
}



void Membrane::calculate_ulm(ArgStruct args){
    int ell_max = args.ell_max;
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
    
    
    if(args.analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    
    
    
    vector<vector< std::complex<double> > > ulm_avg_frame;
    ulm_avg_frame.resize(ell_max+1);
    for (int ell=0; ell<ell_max+1; ell++) {
        ulm_avg_frame[ell].resize(2*ell+1);
        
        for (int m=0; m<2*ell+1; m++) {
            ulm_avg_frame[ell][m].real(0);
            ulm_avg_frame[ell][m].imag(0);
        }
    }
    
    vector<double> membrane_radii_list = get_ulmYlm_vectorlist_for_mesh();
    for (int ell=0; ell<ell_max+1; ell++) {
        //        cout<<ell<<" out of "<<ell_max<<"\r";
        for (int m=-ell; m<ell+1; m++) {
            
            vector<std::complex<double> >  ylm_cc = get_ylm_vectorlist_for_mesh(ell, m, true);
            ulm_avg_frame[ell][m+ell] = calc_vectorlist_vectorlist_surface_integral(ylm_cc, membrane_radii_list);
            
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

void Membrane::calculate_real_ulm(ArgStruct args){
    int ell_max = args.ell_max;
    if(ulm_avg.size() != ell_max+1){
        ulm_avg.clear();
        ulm_std.clear();
        ulm_avg.resize(ell_max+1);
        ulm_std.resize(ell_max+1);
        for (int ell=0; ell<ell_max+1; ell++) {
            ulm_avg[ell].resize(2*ell+1,0);
            ulm_std[ell].resize(2*ell+1,0);
        }
        if (!GenConst::Testmode) {
            cout<<"cleared ulm\n";
        }
        
    }
    
    
    if(args.analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    
    
    vector<vector< double > > ulm_avg_frame;
    ulm_avg_frame.resize(ell_max+1);
    for (int ell=0; ell<ell_max+1; ell++) {
        ulm_avg_frame[ell].resize(2*ell+1);
        
        for (int m=0; m<2*ell+1; m++) {
            ulm_avg_frame[ell][m]=0;
        }
    }
    
    vector<double> membrane_radii_list = get_ulmYlm_vectorlist_for_mesh();
    for (int ell=0; ell<ell_max+1; ell++) {
        //        cout<<ell<<" out of "<<ell_max<<"\r";
        for (int m=-ell; m<ell+1; m++) {
            
            vector<double>  Realylm = get_real_ylm_vectorlist_for_mesh(ell, m);
            ulm_avg_frame[ell][m+ell] = calc_vectorlist_vectorlist_surface_integral(Realylm, membrane_radii_list);
            
        }
    }
    
    
    for (int ell=0; ell<ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            double ulm = ulm_avg_frame[ell][m+ell];
            ulm_avg[ell][m+ell] += ulm*ulm;
            ulm_std[ell][m+ell] += ulm*ulm*ulm*ulm;
            
        }
    }
    
}

void Membrane::write_ulm(ArgStruct args, int file_index){
    std::ofstream wdata;
    wdata.open(args.output_filename[file_index].c_str(), std::ios::app);

    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            
            wdata<<ulm_avg[ell][m+ell]<<"\t";
        }
        wdata<<"\n";
        
    }
    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            
            wdata<<ulm_std[ell][m+ell]<<"\t";
        }
        wdata<<"\n";
    }
//    cout<<args.output_filename<<" written"<<endl;
}
