#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>


using namespace std::complex_literals;



void Membrane::surface_integral_test(){
    
    double theta = M_PI/5.7;
    double phi   = M_PI/3.5;
    int    ell   = 4;
    int    m     = 3;
    std::complex<double> ylm, ylm2, integral;
    std::complex<double> ylm_cc, mult;
    
    
    integral.real(0);
    integral.imag(0);
    for(int i=0;i<Num_of_Nodes;i++){
        ylm = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
        ylm_cc = ylm;
        ylm_cc.imag(-1*imag(ylm));
        
        ylm2 = boost::math::spherical_harmonic(ell+1,m-1,spherical_positions[i][1],spherical_positions[i][2]);
        
        integral+=  ylm_cc* ylm2 * node_voronoi_area[i];
    }
    
    cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell+1<<","<<m-1<<"} = "<<real(integral)*4*M_PI/surface_area_voronoi<<" "<<imag(integral)<<endl;
    
    ell=5;
    m=-2;
    integral.real(0);
    integral.imag(0);
    for(int i=0;i<Num_of_Nodes;i++){
        ylm = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
        ylm_cc = ylm;
        ylm_cc.imag(-1*imag(ylm));
        
        ylm2 = boost::math::spherical_harmonic(ell+1,m-1,spherical_positions[i][1],spherical_positions[i][2]);
        
        integral+=  ylm_cc* ylm2 * node_voronoi_area[i];
    }
    
    cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell+1<<","<<m-1<<"} = "<<real(integral)*4*M_PI/surface_area_voronoi<<" "<<imag(integral)<<endl;
    
    ell=1;
    m=0;
    integral.real(0);
    integral.imag(0);
    for(int i=0;i<Num_of_Nodes;i++){
        ylm = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
        ylm_cc = ylm;
        ylm_cc.imag(-1*imag(ylm));
        
        ylm2 = boost::math::spherical_harmonic(ell+1,m-1,spherical_positions[i][1],spherical_positions[i][2]);
        
        integral+=  ylm_cc* ylm2 * node_voronoi_area[i];
    }
    
    cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell+1<<","<<m-1<<"} = "<<real(integral)*4*M_PI/surface_area_voronoi<<" "<<imag(integral)<<endl;
    
}

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
    
    double voronoi_to_omega_multiplyer = 1./(radius*radius*radius);
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
                double multi = f_theta_phi*voronoi_to_omega_multiplyer*node_voronoi_area[i];
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

void Membrane::calculate_ulm_radiustest(int ell_max, int analysis_averaging_option){
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
    
//    calculate_volume_and_surface_area();
//    double radius = cbrt( 3*volume/(M_PI*4) );
    
//    double radius=0;
//    for (int i=0; i<Num_of_Nodes; i++) {
//        radius += spherical_positions[i][0];
//    }
//    radius/=Num_of_Nodes;
    
    std::complex<double> ylm;
    std::complex<double> ylm_cc;
//    double voronoi_to_omega_multiplyer =1./(radius*radius*radius);
    double voronoi_to_omega_multiplyer =1;
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
                voronoi_to_omega_multiplyer  = 1./(radius*spherical_positions[i][0]*spherical_positions[i][2]);
                double multi = f_theta_phi*voronoi_to_omega_multiplyer*node_voronoi_area[i];
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

void Membrane::calculate_ulm_sub_particles(int ell_max, int analysis_averaging_option){
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
    
    
    std::complex<double> ylm;
    std::complex<double> ylm_cc;
    
    double voronoi_to_omega_multiplyer = 1./(radius*radius*radius);
    
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
    
    vector<int> subparticles {    0,   3,  10,  13,  16,  19,  27,  30,  33,  36,  39,  44,  45,  49,
                                 53,  57,  60,  63,  65,  70,  72,  74,  81,  83,  85,  87,  95,  97,
                                 99, 101, 112, 115, 117, 122, 124, 126, 133, 135, 137, 139, 147, 149,
                                151, 153, 164, 167, 169, 174, 176, 178, 185, 187, 189, 191, 199, 201,
                                203, 205, 216, 219, 221, 226, 228, 230, 237, 239, 241, 243, 251, 253,
                                255, 257, 268, 273, 275, 282, 284, 286, 295, 297, 299, 301, 304, 308,
                                312, 316, 320, 325, 327, 334, 336, 338, 347, 349, 351, 353, 356, 360,
                                364, 368, 372, 377, 379, 386, 388, 390, 399, 401, 403, 405, 408, 412,
                                416, 420, 424, 429, 431, 438, 440, 442, 451, 453, 455, 457, 460, 464,
                                468, 472, 476, 481, 483, 490, 492, 494, 503, 505, 507, 509, 512, 516,
                                520, 524, 528, 533, 535, 542, 544, 546, 555, 557, 559, 561, 572, 577,
                                579, 586, 588, 590, 599, 601, 603, 605, 616, 621, 623, 630, 632, 634,
                                643, 645, 647, 649, 660, 665, 667, 674, 676, 678, 687, 689, 691, 693,
                                704, 709, 711, 718, 720, 722, 731, 733, 735, 737, 748, 753, 755, 762,
                                764, 766, 775, 777, 779, 781, 784, 788, 792, 796, 800, 805, 807, 814,
                                816, 818, 827, 829, 831, 833, 844, 849, 851, 858, 860, 862, 871, 873,
                                875, 877, 888, 893, 895, 902, 904, 906, 915, 917, 919, 921, 932, 937,
                                939, 946, 948, 950, 959, 961, 963, 965, 975, 980, 982, 989, 991, 993};
    
    
    voronoi_to_omega_multiplyer  *= double(Num_of_Nodes)/subparticles.size();
    
    for (int ell=0; ell<ell_max+1; ell++) {
        //        cout<<ell<<" out of "<<ell_max<<"\r";
        for (int m=-ell; m<ell+1; m++) {
            for (int i : subparticles) {
                ylm = boost::math::spherical_harmonic(ell,
                                                      m,
                                                      spherical_positions[i][1],
                                                      spherical_positions[i][2]);
                ylm_cc = ylm;
                ylm_cc.imag(-1*imag(ylm));

                f_theta_phi = (spherical_positions[i][0]-radius);
                double multi = f_theta_phi*voronoi_to_omega_multiplyer*node_voronoi_area[i];
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
    
//    for (int ell=0; ell<ell_max+1; ell++) {
//        for (int m=-ell; m<ell+1; m++) {
//            
//            ulm_avg[ell][m+ell] /= num_frames;
//            
//            ulm_std[ell][m+ell] /= num_frames;
//            ulm_std[ell][m+ell] -= ulm_avg[ell][m+ell]*ulm_avg[ell][m+ell];
//            ulm_std[ell][m+ell] = sqrt(ulm_std[ell][m+ell]/(num_frames-1 ));
//        }
//        
//    }
    
    
    
    
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
}
