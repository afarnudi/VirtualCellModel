//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "Global_functions.hpp"

void Membrane::initialise(std::string Mesh_file_name){
//    T_Kinetic_Energy.resize(100);
    if (!generalParameters.Testmode) {
        cout<<TMEM<<"\nInitialising the Membrane Class..."<<TRESET<<endl;
    }
    
    read_mesh_file(Mesh_file_name);
    
    output_file_neme=Mesh_file_name;
    if (!generalParameters.Testmode) {
        
        cout<<"Number of ... :\n";
        cout<<"Nodes\t"<<Num_of_Nodes<<endl;
        cout<<"Triangles\t"<<Num_of_Triangles<<endl;
    }
    node_mass.resize(Num_of_Nodes, node_global_mass);
    
    update_COM_position();
    Normal_direction_Identifier();
//    Triangle_pair_counter();
    Triangle_pair_identifier();
    if (!generalParameters.Testmode) {
        cout<<"Triangle pairs\t"<<Num_of_Triangle_Pairs<<endl;
    }
    
    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
        cout<<TWARN<<"Warning"<<TRESET<<"! some triangles have less or more neighbour than 3"<<endl;
    }
    Node_Bonds_identifier();
    if (!generalParameters.Testmode) {
        cout<<"Node pairs (Bonds) "<<Num_of_Node_Pairs<<endl;
    }
    
    cout<<"\nOther properties:\n";
    update_average_Membrane_radius();
    cout<<"\nRadius\t"<<Radius<<endl;
//    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
    
    if (AddSphericalHarmonicsMode) {
        generate_undulations();
    }
    if (AddRandomModes) {
        randomundulationgenerator();
    }
    check();
    set_bond_nominal_length();
    set_dihedral_atoms();
    set_bending_nominal_angle();
    set_node_radius();
    calculate_volume_and_surface_area();
    assign_surface_volume_constraints();
    
    if (UseMeanCurvature) {
        mean_curvature_init();
//        export_mesh_properties(Mesh_file_name);
    }
    if (freezeSubLattice) {
        freeze_sub_lattice();
    }
    
    
    if (New_Radius!=-1) {
        check_radius_update_values();
    }
    
    shift_velocity(Shift_velocities_xyzVector[0], Shift_velocities_xyzVector[1], Shift_velocities_xyzVector[2]);
    
    //limiting the labels to 4 charachters for use in the pdb writer
    while (label.length()>3) {
        label.pop_back();
    }
    while (label.length()<3) {
        label += "0";
    }
    if (index>=10){
        label.pop_back();
        label += std::to_string(index);
    } else {
        label += std::to_string(index);
    }
    
//    if (spring_model == potentialModelIndex.Model["FENE"]) {
//        if (FENE_k == 0 || FENE_epsilon == 0 || FENE_max == 0 ) {
//            string errorMessage = TWARN;
//            errorMessage+="Warning";
//            errorMessage+= TRESET;
//            errorMessage+= ". Membrane spring model set to FENE but FENE parameters not set in the membrane configuration file. Please make sure you have set the following parameters: \nFENE_eps\nFENE_k\nFENE_min\nFENE_max (cannot be zero)\n";
//            throw std::runtime_error(errorMessage);
//        }
//    } else
        if (spring_model == potentialModelIndex.Model["None"]) {
        cout<<TWARN<<"\nMembrnae spring model is set to 'None'."<<TRESET<<endl;
    }
    
    if (initial_random_rotation_coordinates){
        if (!generalParameters.Testmode) {
            cout<<TWARN<<"\nRandomly rotating the mesh\n"<<TRESET;
        }
        
        srand (time(NULL));
        double scratch = rand();
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
//        cout<<"theta: "<<theta<<"\n";
//        cout<<"phi  : "<<phi<<"\n";
        rotate_coordinates(theta, phi);
    }
//    if (GenConst::Wantvoronoi){
//        node_voronoi_area.resize(Num_of_Nodes,0);
//    }
    shift_position(Shift_position_xyzVector[0], Shift_position_xyzVector[1], Shift_position_xyzVector[2]);
    if (!generalParameters.Testmode) {
        cout<<TSUCCESS<<"\nMembrane class initiated."<<TRESET<<
                        "\n*************************\n"<<endl;
    }
//    analyse_curvature();
}

void Membrane::analyse_curvature(void){
    cout<<"enter path to xyz:"<<endl;
    string path_to_xyz;
    cin>>path_to_xyz;
    check_if_file_exists(path_to_xyz);
    import_xyz_frames(path_to_xyz);
    load_analysis_coord_frame(0);
//    load_analysis_coord_frame(analysis_coord_frames.size()-1);
    calc_Julicher_Itzykson_bending_props();
    write_vertex_bending_props();
    
    
    double I=0;
    double IB=0;
    double J=0;
    double JV=0;
    double sum_V=0;
    double sum_B=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        I+=Itzykson_numerator[i]/node_voronoi_area[i];
        JV+=Julicher_numerator[i]/node_voronoi_area[i];
        J+=Julicher_numerator[i]/Barycentric_area[i];
        IB+=Itzykson_numerator[i]/Barycentric_area[i];
        sum_V+=node_voronoi_area[i];
        sum_B+=Barycentric_area[i];
    }
    cout<<"I "<<I<<endl;
    cout<<"IB "<<IB<<endl;
    cout<<"J "<<J<<endl;
    cout<<"JV "<<JV<<endl;
    cout<<"V "<<sum_V<<endl;
    cout<<"B "<<sum_B<<endl;
    exit(0);
}

void Membrane::write_vertex_bending_props(){
    
}


void Membrane::calc_Julicher_Itzykson_bending_props(void){
    calculate_surface_area_with_voronoi();
    
    Julicher_numerator.clear();
    Julicher_numerator.resize(Num_of_Nodes,0);
    Itzykson_numerator.clear();
    Itzykson_numerator.resize(Num_of_Nodes,0);
    Barycentric_area.clear();
    Barycentric_area.resize(Num_of_Nodes,0);
    
    for (int node_degree=0; node_degree<nodeOrder_NodeIndex_NodeNeighbourList.size(); node_degree++) {
        for (int node=0; node<nodeOrder_NodeIndex_NodeNeighbourList[node_degree].size(); node++) {
            
            int mem_vertex = nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][0];
            int mem_vertex_j_m1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node].back();
            int mem_vertex_j=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][1];
            int mem_vertex_j_p1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][2];
            
            double phi_ij = M_PI - abs(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
            double ell_ij = vector_length(vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j]));
            phi_ij*=sign_function(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
            Julicher_numerator[mem_vertex]+=ell_ij*phi_ij;
            
            vector<double> r_ij = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex]);
            vector<double> r_ij_p1 = vector_subtract(Node_Position[mem_vertex_j_p1], Node_Position[mem_vertex]);
            Barycentric_area[mem_vertex]+=vector_length(crossvector(r_ij, r_ij_p1));
            
            double cot_sum=0;
            vector<double> sum_Itzykson_numerator_vector(3,0);
            vector<double> r_j_m1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_m1]);
            vector<double> r_j_m1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_m1]);
            cot_sum += innerproduct(r_j_m1_i, r_j_m1_j)/vector_length(crossvector(r_j_m1_i, r_j_m1_j));
            vector<double> r_j_p1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_p1]);
            vector<double> r_j_p1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]);
            cot_sum += innerproduct(r_j_p1_i, r_j_p1_j)/vector_length(crossvector(r_j_p1_i, r_j_p1_j));
            for (int i=0; i<3; i++) {
                sum_Itzykson_numerator_vector[i] += -r_ij[i]*cot_sum;
            }
            
            for (int mem_vertex_neighbour_ind=2; mem_vertex_neighbour_ind<nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node].size()-1; mem_vertex_neighbour_ind++) {
                mem_vertex_j_m1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][mem_vertex_neighbour_ind-1];
                mem_vertex_j=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][mem_vertex_neighbour_ind];
                mem_vertex_j_p1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][mem_vertex_neighbour_ind+1];
                
                phi_ij = M_PI - abs(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
                ell_ij = vector_length(vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j]));
                phi_ij*=sign_function(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
                Julicher_numerator[mem_vertex]+=ell_ij*phi_ij;
                
                r_ij = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex]);
                r_ij_p1 = vector_subtract(Node_Position[mem_vertex_j_p1], Node_Position[mem_vertex]);
                Barycentric_area[mem_vertex]+=vector_length(crossvector(r_ij, r_ij_p1));
                
                cot_sum=0;
                r_j_m1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_m1]);
                r_j_m1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_m1]);
                cot_sum += innerproduct(r_j_m1_i, r_j_m1_j)/vector_length(crossvector(r_j_m1_i, r_j_m1_j));
                r_j_p1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_p1]);
                r_j_p1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]);
                cot_sum += innerproduct(r_j_p1_i, r_j_p1_j)/vector_length(crossvector(r_j_p1_i, r_j_p1_j));
                for (int i=0; i<3; i++) {
                    sum_Itzykson_numerator_vector[i] += -r_ij[i]*cot_sum;
                }
                
            }
            
            mem_vertex_j_m1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node].size()-2];
            mem_vertex_j=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node].back();
            mem_vertex_j_p1=nodeOrder_NodeIndex_NodeNeighbourList[node_degree][node][1];
            
            phi_ij = M_PI - abs(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
            ell_ij = vector_length(vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j]));
            phi_ij*=sign_function(calc_dihedral_angle(Node_Position[mem_vertex_j_m1], Node_Position[mem_vertex_j], Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]));
            Julicher_numerator[mem_vertex]+=ell_ij*phi_ij;
            
            r_ij = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex]);
            r_ij_p1 = vector_subtract(Node_Position[mem_vertex_j_p1], Node_Position[mem_vertex]);
            Barycentric_area[mem_vertex]+=vector_length(crossvector(r_ij, r_ij_p1));
            
            cot_sum=0;
            r_j_m1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_m1]);
            r_j_m1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_m1]);
            cot_sum += innerproduct(r_j_m1_i, r_j_m1_j)/vector_length(crossvector(r_j_m1_i, r_j_m1_j));
            r_j_p1_j = vector_subtract(Node_Position[mem_vertex_j], Node_Position[mem_vertex_j_p1]);
            r_j_p1_i = vector_subtract(Node_Position[mem_vertex], Node_Position[mem_vertex_j_p1]);
            cot_sum += innerproduct(r_j_p1_i, r_j_p1_j)/vector_length(crossvector(r_j_p1_i, r_j_p1_j));
            for (int i=0; i<3; i++) {
                sum_Itzykson_numerator_vector[i] += -r_ij[i]*cot_sum;
            }
            
            
            Julicher_numerator[mem_vertex]*=Julicher_numerator[mem_vertex]*0.125;
            Barycentric_area[mem_vertex]*=1./6.;
            Itzykson_numerator[mem_vertex]=vector_length_squared(sum_Itzykson_numerator_vector)*0.125;
            
            
        }
    }
    
    
}


double Membrane::calc_dihedral_angle(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4){
    vector<double> b1(3,0);
    vector<double> b2(3,0);
    vector<double> b3(3,0);

    double phi=0;
    
    for (int index=0; index<3; index++) {
        b1[index]  = p2[index]-p1[index];
        b2[index]  = p3[index]-p2[index];
        b3[index]  = p4[index]-p3[index];
    }
    
    vector<double> N1 = crossvector(b1, b2);
    vector<double> N2 = crossvector(b2, b3);
    vector<double> b2Normalised = normalise_vector(b2);
    
    vector<double> m1 = crossvector(N1, b2Normalised);
        
    double x = innerproduct(N1, N2);
    double y = innerproduct(m1, N2);
    phi = atan2(y, x);
    return phi;
}

void Membrane::load_analysis_coord_frame(int frame){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] = analysis_coord_frames[frame][i][j];
        }
    }
    
    update_COM_position();
    set_com_to_zero();
    
    update_spherical_positions();
    calculate_dOmega();
    calculate_surface_area_with_voronoi();
}

void Membrane::import_xyz_frames(string path){
    int xyz_Num_of_Nodes = get_xyz_num_of_atoms(path);
    
    std::ifstream read_xyz;
    read_xyz.open(path.c_str());
    
    string line;
    read_xyz.seekg(std::ios::beg);
    
    int num_of_frames = get_xyz_num_of_frames(path, xyz_Num_of_Nodes);
    
    analysis_coord_frames_time.clear();
    analysis_coord_frames_time.resize(num_of_frames);

    analysis_coord_frames.clear();
    analysis_coord_frames.resize(num_of_frames);
    for (int frame_index=0; frame_index<num_of_frames; frame_index++) {
        analysis_coord_frames[frame_index].resize(xyz_Num_of_Nodes);
        
        getline(read_xyz, line);
        getline(read_xyz, line);
        string time=to_string(frame_index);
        if (line.size()!=0){
            auto split = split_and_check_for_comments(line, "Import XYZ: ");
            time = split[1];
        }

        try {
            analysis_coord_frames_time[frame_index]=stod(time);
        } catch (...) {
            analysis_coord_frames_time[frame_index]=frame_index;
        }



        int local_atom_count=0;
        double readCoord[3];
        string readLable;
        
        for (int j=0; j<xyz_Num_of_Nodes; j++) {
            analysis_coord_frames[frame_index][j].resize(3,0);
            read_xyz>>readLable;
            read_xyz>>readCoord[0]>>readCoord[1]>>readCoord[2];
            if (readLable==label) {
                analysis_coord_frames[frame_index][local_atom_count][0]=readCoord[0];
                analysis_coord_frames[frame_index][local_atom_count][1]=readCoord[1];
                analysis_coord_frames[frame_index][local_atom_count][2]=readCoord[2];
                local_atom_count++;
//                cout<<readCoord[0]<<" "<<readCoord[1]<<" "<<readCoord[2]<<endl;
            }

        }
        getline(read_xyz, line);
    }
    
    
}


void Membrane::analysis_init(std::string Mesh_path){
    string mesh_extension = Mesh_path;
    mesh_extension.erase(mesh_extension.begin(),mesh_extension.end()-3 );
    
    if (mesh_extension=="ply"){
        read_ply_file(Mesh_path);
    } else {
        read_gmesh_file(Mesh_path);
    }
    update_average_Membrane_radius();
    
    Normal_direction_Identifier();
    Triangle_pair_counter();
    
    Node_Bonds_identifier();
    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
}


void Membrane::freeze_sub_lattice(void){
    node_mass.resize(0);
    node_mass.resize(Num_of_Nodes,-1);
    for (int i=0; i<nodeOrder_NodeIndex_NodeNeighbourList.size(); i++) {
        for (int j=0; j<nodeOrder_NodeIndex_NodeNeighbourList[i].size(); j++) {
            if (node_mass[nodeOrder_NodeIndex_NodeNeighbourList[i][j][0]]==-1) {
                node_mass[nodeOrder_NodeIndex_NodeNeighbourList[i][j][0]]=node_global_mass;
                for (int k=1; k<nodeOrder_NodeIndex_NodeNeighbourList[i][j].size(); k++) {
                    node_mass[nodeOrder_NodeIndex_NodeNeighbourList[i][j][k]]=0;
                }
            }
            
        }
    }
//    int count_zero=0;
//    int count_mass=0;
//    int count_N=0;
//    
//    for (int i=0; i<Num_of_Nodes; i++) {
//        if (node_mass[i]==0) {
//            count_zero++;
//        } else if (node_mass[i]==node_global_mass){
//            count_mass++;
//        } else if (node_mass[i]==-1){
//            count_N++;
//        }
//    }
//    cout<<"zero "<<count_zero<<endl<<"mass "<<count_mass<<endl<<"negative "<<count_N<<endl;exit(0);
}
