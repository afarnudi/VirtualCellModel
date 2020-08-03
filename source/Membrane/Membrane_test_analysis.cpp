#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>


using namespace std::complex_literals;
using namespace std;


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
    // 252 Subparticles of the 1002 mesh that form a shpere.
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



void Membrane::myWritePDBFrame(int frameNum,
                               std::string traj_name)
{
    
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%6.6f potential energy=%.3f KJ/mole\n",
            0.,
            0.,
            0.);
    int index=0;
    for (int n=0; n< Num_of_Nodes; ++n){
        
        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                "mem0",
                'A',
                double(index),
                Node_Position[n][0],
                Node_Position[n][1],
                Node_Position[n][2],
                0.,
                0.);//,
        //                atoms[n].symbol);
    }
    
    
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}
