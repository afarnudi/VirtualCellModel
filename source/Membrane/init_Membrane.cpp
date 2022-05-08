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
    
    if (UseMeanCurvature) {
        mean_curvature_init();
//        export_mesh_properties(Mesh_file_name);
    }
    
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
    if (CalculateBending) {
        analyse_curvature();
    }
//    loop_ulm_gen();
    
//    calculate_mesh_energy_landscape();
    
}
#include <boost/filesystem.hpp>


void Membrane::loop_ulm_gen(void){
    int L=5;
    int M=3;
    
    string seed = Mesh_file_name;
//    vector<string> results;
//    split(seed, "/", results);
//    seed = results.back();
//    split(seed, "s", results);
//    seed = results.back();
//    split(seed, ".", results);
//    seed=results[0];
    seed="0";
    
    for (double u=0; u<1.00001; u+=0.01) {
//        string propDir = "Results/ULMCurvature/Ordered/amp_"+to_string(u)+"/Lnum_"+to_string(L)+"/Mnum_"+to_string(M)+"/";
        string propDir = "Results/ULMCurvature/Ordered/amp_"+to_string(u)+"/Lnum_"+to_string(L)+"/Mnum_"+to_string(M)+"/"+seed+"/";
        boost::filesystem::create_directories(propDir);
        string propFileName = propDir+"bendingProps.txt";
//        cout<<propFileName<<endl;
        update_spherical_positions();
        for (int i=0; i<Num_of_Nodes; i++) {
            spherical_positions[i][0]=1;
        }
        convert_spherical_positions_to_cartisian();
        add_ulm_mode_real(L, M, u, 1);
        calc_Julicher_Itzykson_bending_props();
        std::ofstream write;
        write.open(propFileName.c_str());
        double I=0;
        double IB=0;
        double J=0;
        double JV=0;
        for (int i=0; i<Num_of_Nodes; i++) {
            I+=Itzykson_numerator[i]/node_voronoi_area[i];
            JV+=Julicher_numerator[i]/node_voronoi_area[i];
            J+=Julicher_numerator[i]/Barycentric_area[i];
            IB+=Itzykson_numerator[i]/Barycentric_area[i];
        }
//        cout<<"Itzykson\t"<<I<<"\n";
        write<<"Itzykson\t"<<I<<"\n";
        write<<"Itzykson-Barycentric\t"<<IB<<"\n";
        write<<"Julicher\t"<<J<<"\n";
        write<<"Julicher-Voronoi\t"<<JV<<"\n";
        
        double curvature = 8*M_PI + u*u*0.5*(L+2)*(L+1)*L*(L-1);
        write<<"Theoretical\t"<<curvature<<"\n";
        write<<"#Itzykson_Numerator Julicher_numerator Voronoi_area Barycentric_area"<<"\n";
        for (int i=0; i<Num_of_Nodes; i++) {
            write<<Itzykson_numerator[i]<<" "<<Julicher_numerator[i]<<" "<<node_voronoi_area[i]<<" "<<Barycentric_area[i]<<"\n";
        }
        write.close();
        write_xyz  (propDir+"coords.xyz");
        myWritePSF(propDir+"coords.psf");
    }
    
    
    exit(0);
    
}
void Membrane::myWritePSF(string traj_name)
{
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile," PSF\n\n");
    fprintf(pFile,"       1  !TITLE\n");
    fprintf(pFile," vmd files (psf,pdb) for VCM\n\n");
    fprintf(pFile,"   %5d !NATOM\n",Num_of_Nodes);
    
    
    for (int n=0; n<Num_of_Nodes; ++n){
        //        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
        fprintf(pFile,"   %5d POLY%5d POLY %4s CEL1 %10.2f    %10.2f 1\n",
                n+1,
                0,
                "Mem0",
                1.0,
                1);
    }
    int num_of_bonds=Num_of_Node_Pairs;
    
    fprintf(pFile,"\n");
    fprintf(pFile,"   %5d !NBOND\n",num_of_bonds);
    int endline_counter=0;
    for (int n=0; n<Num_of_Node_Pairs; ++n){
//        if (bonds[n].type != potentialModelIndex.Model["None"]) {
            fprintf(pFile,"   %5d   %5d",
                    Node_Bond_list[n][0]+1,
                    Node_Bond_list[n][1]+1);
            endline_counter++;
            if (endline_counter>3) {
                endline_counter=0;
                fprintf(pFile,"\n");
            }
//        }
        
    }
    fclose (pFile);
}
void Membrane::write_xyz  (string traj_name)
{
    ofstream writexyz(traj_name.c_str());
    writexyz<<Num_of_Nodes<<"\n";
    writexyz<<"timePs ";
    writexyz<<0<<setprecision(16);
    writexyz<<" potential_energy_inKJpermol ";
    writexyz<<0<<setprecision(16);
    writexyz<<" energy_inKJpermol ";
    writexyz<<0<<setprecision(16)<<"\n";
    for (int n=0; n<Num_of_Nodes; n++) {
        writexyz<<"Mem0\t"<<Node_Position[n][0]<<"\t"<<Node_Position[n][1]<<"\t"<<Node_Position[n][2]<<setprecision(16)<<"\n";
    }
    writexyz.close();
}

void Membrane::analyse_curvature(void){
    
    bool simulationBeginning=true;
    if (xyzPostAnalysisPath!="Path/to/my/xyzfile.xyz") {
        simulationBeginning=false;
        check_if_file_exists(xyzPostAnalysisPath);
        import_xyz_frames(xyzPostAnalysisPath);
        if (xyzPostAnalysisPathFrame==-1) {
            xyzPostAnalysisPathFrame =analysis_coord_frames.size()-1;
        }
        load_analysis_coord_frame(xyzPostAnalysisPathFrame);
        calc_Julicher_Itzykson_bending_props();
    }
    
    calc_Julicher_Itzykson_bending_props();
    
    
    
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
    cout<<"Itzykson              "<<I<<endl;
    cout<<"Itzykson-Barycentric  "<<IB<<endl;
    cout<<"Julicher              "<<J<<endl;
    cout<<"Julicher-Voronoi      "<<JV<<endl;
    cout<<"\nVoronoi area          "<<sum_V<<endl;
    cout<<"Barycentric area      "<<sum_B<<endl;
    if (!simulationBeginning) {
        write_vertex_bending_props();
        exit(0);
    } else {
        xyzPostAnalysisPathFrame=0;
    }
    
}

void Membrane::write_vertex_bending_props(){
    string propFileName = generalParameters.trajectory_file_name+"_"+label+"_BendingProps_frame_"+to_string(xyzPostAnalysisPathFrame)+".txt";
//    cout<<propFileName<<endl;exit(0);
    std::ofstream write;
    write.open(propFileName.c_str(), std::ios::app);
    double I=0;
    double IB=0;
    double J=0;
    double JV=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        I+=Itzykson_numerator[i]/node_voronoi_area[i];
        JV+=Julicher_numerator[i]/node_voronoi_area[i];
        J+=Julicher_numerator[i]/Barycentric_area[i];
        IB+=Itzykson_numerator[i]/Barycentric_area[i];
    }
    write<<"Itzykson\t"<<I<<endl;
    write<<"Itzykson-Barycentric\t"<<IB<<endl;
    write<<"Julicher\t"<<J<<endl;
    write<<"Julicher-Voronoi\t"<<JV<<endl;
    write<<"#Itzykson_Numerator Julicher_numerator Voronoi_area Barycentric_area"<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write<<Itzykson_numerator[i]<<" "<<Julicher_numerator[i]<<" "<<node_voronoi_area[i]<<" "<<Barycentric_area[i]<<endl;
    }
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





/**Energy Landscape calculations:
 This code will produce data required by a post analysis Python script.
 It calculates the energy change due to translational movment of one node on the surface of a sphere.
*/

void Membrane::calculate_mesh_energy_landscape(void){
//    int node_index=14;
    vector<vector<double> > node_pos_copy = Node_Position;
    
    
    for (int node_index=13; node_index<14; node_index++) {
        //Map on unit sphere
        Node_Position = node_pos_copy;
        update_spherical_positions();
        for (int i=0; i<Num_of_Nodes; i++) {
            spherical_positions[i][0]=1;
        }
        convert_spherical_positions_to_cartisian();
        
        //set node to north pole
        vector<double> r_theta_phi;
        r_theta_phi = convert_cartesian_to_spherical(Node_Position[node_index][0],
                                                     Node_Position[node_index][1],
                                                     Node_Position[node_index][2]);
        double theta0 = r_theta_phi[1];
        double phi0   = r_theta_phi[2];
        rotate_coordinates(-theta0, -phi0);
        
        double range=0.3;
        int segments=40;
        
        vector<vector<double> > J_curvature(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > I_curvature(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > JV_curvature(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > IB_curvature(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > B_area(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > V_area(2*segments, vector<double>(2*segments,0));
        vector<vector<double> > temp_volume(2*segments, vector<double>(2*segments,0));
    
        cout<<node_index<<" "<<flush;
        for (int j=0; j<2*segments; j++) {
            cout<<j<<" "<<flush;
            Node_Position[node_index][1]=(1-j/double(segments))*range;
            for (int i=-segments; i<segments; i++) {
                Node_Position[node_index][0]=i*range/segments;
                Node_Position[node_index][2]=sqrt(1-Node_Position[node_index][0]*Node_Position[node_index][0]-Node_Position[node_index][1]*Node_Position[node_index][1]);
    //            cout<<"3"<<endl;
                
//                analyse_curvature();
                calculate_volume_and_surface_area();
                temp_volume[2*segments-j-1][i+segments]+=volume;
                calc_Julicher_Itzykson_bending_props();
    //            cout<<"4"<<endl;
                for (int nid=0; nid<Num_of_Nodes; nid++) {
                    J_curvature[2*segments-j-1][i+segments]+=Julicher_numerator[nid]/Barycentric_area[nid];
                    I_curvature[2*segments-j-1][i+segments]+=Itzykson_numerator[nid]/node_voronoi_area[nid];
                    IB_curvature[2*segments-j-1][i+segments]+=Itzykson_numerator[nid]/Barycentric_area[nid];
                    JV_curvature[2*segments-j-1][i+segments]+=Julicher_numerator[nid]/node_voronoi_area[nid];
                    B_area[2*segments-j-1][i+segments]+=Barycentric_area[nid];
                    V_area[2*segments-j-1][i+segments]+=node_voronoi_area[nid];
                }

            }
        }
//        cout<<"2"<<endl;
        string JName= "Node_curvatures/Node_"+to_string(node_index)+"_Julicher_curvature.txt";
        string IName= "Node_curvatures/Node_"+to_string(node_index)+"_Itzykson_curvature.txt";
        string IBName= "Node_curvatures/Node_"+to_string(node_index)+"_ItzyksonBarycentric_curvature.txt";
        string JVName= "Node_curvatures/Node_"+to_string(node_index)+"_JulicherVoronoi_curvature.txt";
        string BName= "Node_curvatures/Node_"+to_string(node_index)+"_Barycentric_area.txt";
        string VName= "Node_curvatures/Node_"+to_string(node_index)+"_Voronoi_area.txt";
        string volumeName= "Node_curvatures/Node_"+to_string(node_index)+"_volume.txt";
        ofstream wJ(JName.c_str());
        ofstream wI(IName.c_str());
        ofstream wIB(IBName.c_str());
        ofstream wJV(JVName.c_str());
        ofstream wB(BName.c_str());
        ofstream wV(VName.c_str());
        ofstream wvol(volumeName.c_str());
        for (int i=0; i<2*segments; i++) {
            for (int j=0; j<2*segments; j++) {
                wJ<<setprecision(10)<<J_curvature[i][j]<<" ";
                wI<<setprecision(10)<<I_curvature[i][j]<<" ";
                wIB<<setprecision(10)<<IB_curvature[i][j]<<" ";
                wJV<<setprecision(10)<<JV_curvature[i][j]<<" ";
                wB<<setprecision(10)<<B_area[i][j]<<" ";
                wV<<setprecision(10)<<V_area[i][j]<<" ";
                wvol<<setprecision(10)<<temp_volume[i][j]<<" ";
            }
            wJ<<endl;
            wI<<endl;
            wIB<<endl;
            wJV<<endl;
            wB<<endl;
            wV<<endl;
            wvol<<endl;
        }
    }
    
    
    exit(0);
}


