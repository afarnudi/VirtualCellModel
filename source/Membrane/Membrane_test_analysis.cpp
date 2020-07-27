#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>


using namespace std::complex_literals;



void Membrane::surface_integral_test(){
    
    for (int j =0; j<1000; j++) {
        int    ell   = rand()%20;
        int    m     = rand()%(ell+1)*(rand()%2)*(-1);
        int    ell2   = rand()%20;
        int    m2     = rand()%(ell+1)*(rand()%2)*(-1);
        
        std::complex<double> ylm, ylm2, ylm3, integral, integral_same;
        std::complex<double> ylm_cc, mult;
        
        
        integral.real(0);
        integral.imag(0);
        
        integral_same.real(0);
        integral_same.imag(0);
        for(int i=0;i<Num_of_Nodes;i++){
            ylm = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
            ylm_cc = ylm;
            ylm_cc.imag(-1*imag(ylm));
            
            ylm2 = boost::math::spherical_harmonic(ell+1,m-1,spherical_positions[i][1],spherical_positions[i][2]);
            
            ylm3 = boost::math::spherical_harmonic(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
            
            integral+=  ylm_cc* ylm2 * node_voronoi_area[i];
            integral_same+=  ylm_cc* ylm3 * node_voronoi_area[i];
        }
        
        integral *= 4*M_PI/surface_area_voronoi;
        integral_same *= 4*M_PI/surface_area_voronoi;
        if ( real(integral)>0.01  ){
            if (ell != ell2 && m != m2) {
                cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell2<<","<<m2<<"} = "<<real(integral)<<" Imag: "<<imag(integral)<<endl;
            }
            
        }
        if ( (real(integral_same)-1)>0.01 || (real(integral_same)-1)<-0.01 ){
            
            cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell<<","<<m<<"} = "<<real(integral_same)<<" Imag: "<<imag(integral_same)<<endl;
        }
        
    }
    
    
    
}

void Membrane::surface_integral_test_real(){
    
    for (int j=0; j<1000; j++) {
        int    ell   = rand()%20;
        int    m     = rand()%(ell+1)*(rand()%2)*(-1);
        int    ell2   = rand()%20;
        int    m2     = rand()%(ell+1)*(rand()%2)*(-1);
        
        double ylm, ylm2, ylm3, integral=0, integral_same=0;
        double ylm_cc, mult;
        
        
        for(int i=0;i<Num_of_Nodes;i++){
            ylm = real_harmonics(ell, m, spherical_positions[i][1], spherical_positions[i][2]);
            ylm_cc = ylm;
            
            ylm2 = real_harmonics(ell2, m2, spherical_positions[i][1], spherical_positions[i][2]);
            ylm3 = real_harmonics(ell, m, spherical_positions[i][1], spherical_positions[i][2]);
            integral+=  ylm_cc* ylm2 * node_voronoi_area[i];
            integral_same+=  ylm_cc* ylm3 * node_voronoi_area[i];
        }
        integral *= 4*M_PI/surface_area_voronoi;
        integral_same *= 4*M_PI/surface_area_voronoi;
        if ( integral>0.01  ){
            if (ell != ell2 && m != m2) {
                cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell2<<","<<m2<<"} = "<<integral<<endl;
            }
            
        }
        if ( (integral_same-1)>0.01 || (integral_same-1)<-0.01 ){
            cout<<"Integral y_{"<<ell<<","<<m<<"}* x y_{"<<ell<<","<<m<<"} = "<<integral_same<<endl;
        }
        
    }
    
    
}

using namespace std;


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
    
    
    //    double radius = sqrt( surface_area_voronoi/(M_PI*4) );
    
    //        calculate_volume_and_surface_area();
    //        double radius = cbrt( 3*volume/(M_PI*4) );
    //
    double radius=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        radius += spherical_positions[i][0];
    }
    radius/=Num_of_Nodes;
//    cout<<"radius = "<<radius<<endl<<endl;
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
                voronoi_to_omega_multiplyer  = 1./(radius*spherical_positions[i][0]*spherical_positions[i][0]);
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


void Membrane::generate_ulm_mode(int Ell, int M, double uLM){
    if (Ell<0){
        cout<<"\nl must be positive"<<endl;
        exit(EXIT_FAILURE);
    }
    if (abs(M)>Ell){
        cout<<"\n |m| cannot be larger than l"<<endl;
        exit(EXIT_FAILURE);
    }
    
    std::complex<double> ylm;
    //    std::complex<double> ylm_cc;
    double f_theta_phi;
    
    double radius =10;
    
    for(int i=0;i<Num_of_Nodes;i++){
        ylm = boost::math::spherical_harmonic(Ell,M,spherical_positions[i][1],spherical_positions[i][2]);
        //        ylm_cc = ylm;
        //        ylm_cc.imag(-1*imag(ylm));
        
        f_theta_phi = radius*uLM*real(ylm);
        spherical_positions[i][0] = radius+f_theta_phi;
        //        Node_Position[i]=convert_spherical_to_cartesian(spherical_positions[i][0], spherical_positions[i][1], spherical_positions[i][2]);
        Node_Position[i][0]=spherical_positions[i][0]*sin(spherical_positions[i][1])
        *cos(spherical_positions[i][2]);
        Node_Position[i][1]=spherical_positions[i][0]*sin(spherical_positions[i][1])
        *sin(spherical_positions[i][2]);
        Node_Position[i][2]=spherical_positions[i][0]*cos(spherical_positions[i][1]);
        
    }
    
}

void Membrane::generate_ulm_mode_real(int Ell, int M, double uLM){
    
    if (Ell<0){
        cout<<"\nl must be positive"<<endl;
        exit(EXIT_FAILURE);
    }
    if (abs(M)>Ell){
        cout<<"\n |m| cannot be larger than l"<<endl;
        exit(EXIT_FAILURE);
    }
    
    
    double radius = 1;
    
    for(int i=0;i<Num_of_Nodes;i++){
        double ylm = real_harmonics(Ell, M, spherical_positions[i][1],spherical_positions[i][2]);
        double ylm2 = real_harmonics(8, 4, spherical_positions[i][1],spherical_positions[i][2]);
        double f_theta_phi = radius*uLM*ylm;
        spherical_positions[i][0] = radius + f_theta_phi + radius*uLM*ylm2;
        //        Node_Position[i]=convert_spherical_to_cartesian(spherical_positions[i][0], spherical_positions[i][1], spherical_positions[i][2]);
        Node_Position[i][0]=spherical_positions[i][0]*sin(spherical_positions[i][1])
        *cos(spherical_positions[i][2]);
        Node_Position[i][1]=spherical_positions[i][0]*sin(spherical_positions[i][1])
        *sin(spherical_positions[i][2]);
        Node_Position[i][2]=spherical_positions[i][0]*cos(spherical_positions[i][1]);
        
    }
    
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

double Membrane::real_harmonics(int Ell, int M, double theta, double phi){
    std::complex<double> ylm;
    
    double realYLM=0;
    
    if (M<0){
        ylm = boost::math::spherical_harmonic(Ell,-M,theta,phi);
        
        realYLM = imag(ylm)*pow(-1,M)*sqrt(2);
        
    } else if (M==0){
        ylm = boost::math::spherical_harmonic(Ell,0,theta,phi);
        
        realYLM = real(ylm);
    } else {
        
        ylm = boost::math::spherical_harmonic(Ell,M,theta,phi);
        
        realYLM = real(ylm)*pow(-1,M)*sqrt(2);
        
    }
    return realYLM;
}


void Membrane::calculate_ulm_radiustest_real(int ell_max, int analysis_averaging_option){
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
    
    double radius=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        radius += spherical_positions[i][0];
    }
    radius/=Num_of_Nodes;
//    cout<<"radius = "<<radius<<endl<<endl;
    double ylm;
//    double ylm_cc;
    //    double voronoi_to_omega_multiplyer =1./(radius*radius*radius);
    double voronoi_to_omega_multiplyer =1;
    double f_theta_phi;
    
    
    vector<vector< double > > ulm_avg_frame;
    ulm_avg_frame.resize(ell_max+1);
    for (int ell=0; ell<ell_max+1; ell++) {
        ulm_avg_frame[ell].resize(2*ell+1);
        
        for (int m=0; m<2*ell+1; m++) {
            ulm_avg_frame[ell][m]=0;
        }
    }
    
    for (int ell=0; ell<ell_max+1; ell++) {
        //        cout<<ell<<" out of "<<ell_max<<"\r";
        for (int m=-ell; m<ell+1; m++) {
            for(int i=0;i<Num_of_Nodes;i++){
                ylm = real_harmonics(ell,m,spherical_positions[i][1],spherical_positions[i][2]);
//                ylm_cc = ylm;
                
                f_theta_phi = (spherical_positions[i][0]-radius);
                voronoi_to_omega_multiplyer  = 1./(spherical_positions[i][0]*spherical_positions[i][0]);
                double d_omega = voronoi_to_omega_multiplyer*node_voronoi_area[i];
//                ylm_cc =ylm * multi;
                
                ulm_avg_frame[ell][m+ell] += ylm * f_theta_phi * d_omega/radius;
            }
            
        }
    }
    
    for (int ell=0; ell<ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            double ulm = ulm_avg_frame[ell][m+ell];
            ulm_avg[ell][m+ell] += ulm*ulm;
            ulm_std[ell][m+ell] += ulm*ulm*ulm*ulm;
            
        }
    }
    
    //    cout<<endl;
}
