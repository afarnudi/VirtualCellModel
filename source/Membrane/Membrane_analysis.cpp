#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>

void Membrane::rotate_coordinates(double theta, double phi){
    vector<vector<double> > Rotationy;
    vector<vector<double> > Rotationz;
    
    Rotationy.resize(3);
    Rotationz.resize(3);
    for (int i=0; i<3; i++) {
        Rotationy[i].resize(3);
        Rotationz[i].resize(3);
    }
    if (theta != 0) {
        Rotationy[0][0] = cos(theta);
        Rotationy[0][1] = 0;
        Rotationy[0][2] = sin(theta);
        
        Rotationy[1][0] = 0;
        Rotationy[1][1] = 1;
        Rotationy[1][2] = 0;
        
        Rotationy[2][0] = -sin(theta);
        Rotationy[2][1] = 0;
        Rotationy[2][2] = cos(theta);
    }
    
    if (phi !=0) {
        Rotationz[0][0] = cos(phi);
        Rotationz[0][1] = -sin(phi);
        Rotationz[0][2] = 0;
        
        Rotationz[1][0] = sin(phi);
        Rotationz[1][1] = cos(phi);
        Rotationz[1][2] = 0;
        
        Rotationz[2][0] = 0;
        Rotationz[2][1] = 0;
        Rotationz[2][2] = 1;
    }
    
    
    for (int i=0; i<Num_of_Nodes; i++) {
        double x = Node_Position[i][0];
        double y = Node_Position[i][1];
        double z = Node_Position[i][2];
        
        if (phi != 0) {
            Node_Position[i][0] = x*Rotationz[0][0] + y*Rotationz[0][1] + z*Rotationz[0][2];
            Node_Position[i][1] = x*Rotationz[1][0] + y*Rotationz[1][1] + z*Rotationz[1][2];
            Node_Position[i][2] = x*Rotationz[2][0] + y*Rotationz[2][1] + z*Rotationz[2][2];
        }
        
        if (theta !=0) {
            x = Node_Position[i][0];
            y = Node_Position[i][1];
            z = Node_Position[i][2];
            
            Node_Position[i][0] = x*Rotationy[0][0] + y*Rotationy[0][1] + z*Rotationy[0][2];
            Node_Position[i][1] = x*Rotationy[1][0] + y*Rotationy[1][1] + z*Rotationy[1][2];
            Node_Position[i][2] = x*Rotationy[2][0] + y*Rotationy[2][1] + z*Rotationy[2][2];
        }
        
        
        
        
    }
    
    
    
}



double Membrane::calc_theta_angle_ABC(int node_A, int node_B, int node_C){
    /**A is the is the middle point of the angle*/
    double AB[3], AC[3], ACxAB[3];
    AB[0] = Node_Position[node_B][0] - Node_Position[node_A][0];
    AB[1] = Node_Position[node_B][1] - Node_Position[node_A][1];
    AB[2] = Node_Position[node_B][2] - Node_Position[node_A][2];
    
    AC[0] = Node_Position[node_C][0] - Node_Position[node_A][0];
    AC[1] = Node_Position[node_C][1] - Node_Position[node_A][1];
    AC[2] = Node_Position[node_C][2] - Node_Position[node_A][2];
    
    crossvector(ACxAB, AB, AC);
    return innerproduct(AB, AC)/vector_length(ACxAB);
}



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
#include <boost/geometry.hpp>

using namespace std::complex_literals;
namespace bg = boost::geometry;

vector<double> convert_cartesian_to_spherical(double x, double y, double z){
    vector<double> r_theta_phi;
    r_theta_phi.resize(3,0);
    
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > p1;
    bg::model::point<double, 3, bg::cs::cartesian> p3(x,y,z);
    bg::transform(p3, p1);
    
    r_theta_phi[0] = p1.get<2>();
    r_theta_phi[1] = p1.get<1>();
    r_theta_phi[2] = p1.get<0>();
    return r_theta_phi;
}
void Membrane::calculate_surface_area_with_voronoi(){
    
    int triangle1_node_A,
    triangle1_node_B,
    triangle1_node_C,
    triangle2_node_D;
    
    vector<vector<double> > cot_theta_list;
    cot_theta_list.resize(Num_of_Node_Pairs);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        cot_theta_list[i].resize(2,0);
        
        triangle1_node_C = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][0];
        triangle1_node_B = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][1];
        triangle1_node_A = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][2];
        triangle2_node_D = Triangle_Pair_Nodes[ Bond_triangle_neighbour_indices[i] ][3];
        
        cot_theta_list[i][0] = calc_theta_angle_ABC(triangle1_node_C, triangle1_node_A, triangle1_node_B );
        cot_theta_list[i][1] = calc_theta_angle_ABC(triangle2_node_D, triangle1_node_A, triangle1_node_B );
    }
    
    surface_area_voronoi=0;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        node_voronoi_area[i]=0;
        for (int j=0; j<Node_neighbour_list[i].size(); j++) {
            //            cout<<"Node_neighbour_list["<<i<<"]["<<j<<"]="<<node_pair_vec[i][j]<<endl;
            int node_1 = i;
            int node_2 = Node_neighbour_list[node_1][j];
            int bond12 = Node_neighbour_list_respective_bond_index[node_1][j];
            
            double bond_vec[3]={Node_Position[node_1][0]-Node_Position[ node_2 ][0],
                Node_Position[node_1][1]-Node_Position[ node_2 ][1],
                Node_Position[node_1][2]-Node_Position[ node_2 ][2]};
            node_voronoi_area[i] += 0.125 * vector_length(bond_vec)* vector_length(bond_vec)* (cot_theta_list[bond12 ][0] +
                                                                                               cot_theta_list[bond12 ][1]);
        }
        surface_area_voronoi += node_voronoi_area[i];
        
        
    }
}

void Membrane::update_spherical_positions(){
    spherical_positions.clear();
    spherical_positions.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        spherical_positions[i].resize(3,0);
        spherical_positions[i] = convert_cartesian_to_spherical(Node_Position[i][0],                                                                       Node_Position[i][1],
                                                                Node_Position[i][2]);
    }
}

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


void Membrane::set_com_to_zero(){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] -= COM_position[j];
        }
    }
}

void Membrane::rotate_particle_to_axes(void){
    
    //find the angles of the node that used to be on the z axis
    vector<double> r_theta_phi;
    r_theta_phi = convert_cartesian_to_spherical(Node_Position[z_node_index][0],
                                                 Node_Position[z_node_index][1],
                                                 Node_Position[z_node_index][2]);
    
    double theta0 = r_theta_phi[1];
    double phi0   = r_theta_phi[2];
    
    rotate_coordinates(-theta0, -phi0);
    //rotate the particles so the node that was on the z axis at the beginning of the run will lie on the z axis
    r_theta_phi = convert_cartesian_to_spherical(Node_Position[y_node_index][0],
                                                 Node_Position[y_node_index][1],
                                                 Node_Position[y_node_index][2]);
    
    phi0   = r_theta_phi[2];
    //rotate the particles so the node that was on the y axis at the beginning of the run will lie on the y axis
    rotate_coordinates(0, M_PI-phi0);
}

void Membrane::find_nodes_on_the_z_and_y_axis(){
    for (int i=0; i<Num_of_Nodes; i++) {
        if ( (abs( pdb_frames[0][i][0]-1 )<0.01) &&
             (abs( pdb_frames[0][i][1]-1 )<0.01) ) {
            z_node_index=i;
        }
        if ( (abs( pdb_frames[0][i][0]-rescale_factor )<0.01) &&
             (abs( pdb_frames[0][i][2]-rescale_factor )<0.01) ) {
            y_node_index=i;
        }
        if( (y_node_index != -1) && (z_node_index != -1)){
            cout<<"Founds nodes lying on the z and y axis: "<<endl;
            cout<<"z nodes "<<z_node_index<<"\ny node "<<y_node_index<<endl;
            break;
        }
    }
    if(y_node_index == -1) {
        cout<<"Could not find a node on the y axis"<<endl;
    }
    if( (y_node_index == -1) || (z_node_index == -1)){
        cout<<"Could not find a node on the z axis"<<endl;
    }
}

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
            find_nodes_on_the_z_and_y_axis();
        } else {
            cout<<"using:\n z_node_index "<<z_node_index<<"\ny_node_index"<<y_node_index<<endl;
        }
        rotate_particle_to_axes();
    }
    
    
    update_spherical_positions();
    calculate_surface_area_with_voronoi();
}


int Membrane::import_pdb_frames(string filename){
    
    
    int num_of_frames = count_pdb_frames(filename, Num_of_Nodes);
//    cout<<"I'm out\n";
    ifstream read_pdb;
    read_pdb.open(filename.c_str());
    if (!read_pdb) {
        cout << "Unable to read "<<filename<<endl;
        exit(EXIT_FAILURE);
    }
    string line;
    read_pdb.seekg(std::ios::beg);
    
    cout<<"num of frames = "<<num_of_frames<<endl;
    pdb_frames.resize(num_of_frames);
    for (int i=0; i<num_of_frames; i++) {
        pdb_frames[i].resize(Num_of_Nodes);
        for (int j=0; j<Num_of_Nodes; j++) {
            pdb_frames[i][j].resize(3);
            pdb_frames[i][j][0]=0;
            pdb_frames[i][j][1]=0;
            pdb_frames[i][j][2]=0;
        }
    }

    for (int frame_index= 0; frame_index<num_of_frames; frame_index++){

        getline(read_pdb, line);
        getline(read_pdb, line);

        for (int node_index= 0; node_index<Num_of_Nodes; node_index++) {
            
            getline(read_pdb, line);
            istringstream iss(line);
            vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
            
            if(split.size()==11){
                pdb_frames[frame_index][node_index][0] = stod(split[6]);
                pdb_frames[frame_index][node_index][1] = stod(split[7]);
                pdb_frames[frame_index][node_index][2] = stod(split[8]);
            } else if(split.size()==10){
                if (split[6].size()<=9){
                    
                    pdb_frames[frame_index][node_index][0] = stod(split[6]);
                    int it=-1;
                    for(int i=0; i<split[7].size();i++){
                        if (split[7][i] == '.' ){
                            it=i+4;
                            break;
                        }
                    }
                    
                    
                    string coor;
                    coor = split[7];
                    coor.erase(coor.begin() + it, coor.end());
                    pdb_frames[frame_index][node_index][1] = stod(coor);
                    
                    coor = split[7];
                    coor.erase(coor.begin() + 0,coor.begin()+ it );
                    pdb_frames[frame_index][node_index][2] = stod(coor);
                    
                    
                    
                } else {
                    int it=-1;
                    for(int i=0; i<split[6].size();i++){
                        if (split[6][i] == '.' ){
                            it=i+4;
                            break;
                        }
                    }
                    
                    string coor;
                    coor = split[6];
                    coor.erase(coor.begin() + it, coor.end());
                    pdb_frames[frame_index][node_index][0] = stod(coor);
                    
                    coor = split[6];
                    coor.erase(coor.begin() + 0,coor.begin()+ it );
                    pdb_frames[frame_index][node_index][1] = stod(coor);
                    pdb_frames[frame_index][node_index][2] = stod(split[7]);
                }
            } else if(split.size()==9){
                int it[3] = {-1,-1,-1};
                int it_index=0;
                for(int i=0; i<split[6].size();i++){
                    if (split[6][i] == '.' ){
                        it[it_index]=i+4;
                        it_index++;
                    }
                }
                
                string coor;
                coor = split[6];
                coor.erase(coor.begin() + it[0], coor.end());
                pdb_frames[frame_index][node_index][0] = stod(coor);
                
                coor = split[6];
                coor.erase(coor.begin() + it[1], coor.end());
                coor.erase(coor.begin() + 0,coor.begin()+ it[0] );
                pdb_frames[frame_index][node_index][1] = stod(coor);
                
                coor = split[6];
                coor.erase(coor.begin() + 0,coor.begin()+ it[1] );
                pdb_frames[frame_index][node_index][2] = stod(coor);
                
            }

        }
        getline(read_pdb, line);
    }
    return num_of_frames;
}

void Membrane::pdb_to_bin(string filename){
//    ifstream read_pdb;
//    read_pdb.open(filename.c_str());
//    if (!read_pdb) {
//        cout << "Unable to read "<<filename<<endl;
//        exit(EXIT_FAILURE);
//    }
//    string line;
//    string write_name = filename;
//    write_name.erase(write_name.end()-4,write_name.end());
//    write_name += ".bin";
//    ofstream writebin("student.dat", ios::out | ios::binary);
//    if(!writebin) {
//        cout << "Cannot write file "<<write_name << endl;
//        exit(EXIT_FAILURE);
//    }
//
//
//    read_pdb.seekg(std::ios::beg);
//    while(getline(read_pdb, line)){
//        getline(read_pdb, line);
//        cout<<line<<endl;
//        for (int node_index= 0; node_index<Num_of_Nodes; node_index++) {
//            getline(read_pdb, line);
//            istringstream iss(line);
//            vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
//
//            if(split.size()==11){
//                Node_Position[node_index][0] = stod(split[6]);
//                Node_Position[node_index][1] = stod(split[7]);
//                Node_Position[node_index][2] = stod(split[8]);
//            } else if(split.size()==10){
//                if (split[6].size()<=9){
//
//                    Node_Position[node_index][0] = stod(split[6]);
//                    int it=0;
//                    for(int i=1; i<split[7].size();i++){
//                        if (split[7][i] == '-' ){
//                            it=i;
//                        }
//                    }
//                    string coor;
//                    coor = split[7];
//                    coor.erase(coor.begin() + it, coor.end());
//                    Node_Position[node_index][1] = stod(coor);
//
//                    coor = split[7];
//                    coor.erase(coor.begin() + 0,coor.begin()+ it );
//                    Node_Position[node_index][2] = stod(coor);
//
//                } else {
//                    int it=0;
//                    for(int i=1; i<split[6].size();i++){
//                        if (split[6][i] == '-' ){
//                            it=i;
//                        }
//                    }
//                    string coor;
//                    coor = split[6];
//                    coor.erase(coor.begin() + it, coor.end());
//                    Node_Position[node_index][0] = stod(coor);
//
//                    coor = split[6];
//                    coor.erase(coor.begin() + 0,coor.begin()+ it );
//                    Node_Position[node_index][1] = stod(coor);
//                    Node_Position[node_index][2] = stod(split[7]);
//                }
//            } else if(split.size()==9){
//                int it[2] = {0,0};
//                int it_index=0;
//                for(int i=1; i<split[6].size();i++){
//                    if (split[6][i] == '-' ){
//                        it[it_index]=i;
//                        it_index++;
//                    }
//                }
//                string coor;
//                coor = split[6];
//                coor.erase(coor.begin() + it[0], coor.end());
//                Node_Position[node_index][0] = stod(coor);
//
//                coor = split[6];
//                coor.erase(coor.begin() + it[1], coor.end());
//                coor.erase(coor.begin() + 0,coor.begin()+ it[0] );
//                Node_Position[node_index][1] = stod(coor);
//
//                coor = split[6];
//                coor.erase(coor.begin() + 0,coor.begin()+ it[1] );
//                Node_Position[node_index][2] = stod(coor);
//
//            }
//        }
//        getline(read_pdb, line);
//    }
//    //    for (int i=0; i<split.size(); i++) {
//    //        cout<<split[i]<<endl;
//    //    }
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
            ulm_avg[ell][m+ell] += ulm;
            ulm_std[ell][m+ell] += ulm*ulm;
            
        }
    }
    
    //    cout<<endl;
}

void Membrane::write_ulm(int ell_max, string traj_name, double num_frames){
    
    traj_name.erase(traj_name.end()-4,traj_name.end() );
    string traj_file_name = traj_name + "_ulmt_cpp.txt";
    
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
