//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_constants.h"
#include "General_functions.hpp"

using std::ofstream;


void Membrane:: write_pov_traj(std::string traj_name, std::string label, int currentStep){
    ///=============GENRAL===============
    int w1,w2;
    //////________
    std::string pov_file_name;
    pov_file_name=traj_name;
    pov_file_name=pov_file_name+"_of_step_";
    pov_file_name=pov_file_name+std::to_string((currentStep));
    pov_file_name=pov_file_name+".pov";
    
    ofstream pov;
    pov.open(pov_file_name.c_str() );
    pov << std:: fixed;
    
    /////_________
    pov<< "global_settings {   ambient_light rgb<1, 1, 1>   photons {  spacing 0.2  autostop 1  media 60  max_trace_level 6 } } "<<endl;
    pov<< " #include \"colors.inc\"   " <<endl;
    pov<< " #include \"textures.inc\"   " <<endl;
    pov<< " #include \"metals.inc\"   " <<endl;
    pov<< " #include \"glass.inc\"   " <<endl;
    pov<< " #default{ finish{ ambient 0.1 diffuse 0.9 }}   " <<endl;
    
    pov<< "  camera{  location<50,50,50>  look_at <-4,-20,0>    right 0.5*4/3*x     up 0.5*y}"<<endl;
    pov << "background {White} light_source { <000, 100, 000> color White }   light_source { <100, 100, 100> color White }   light_source { <-100, 100, 100> color White } light_source { <00, 100, 100> color White } light_source { <00, 200, 00> color White }"<<endl;
    
    
    
    pov<< "#declare membrane = texture {    pigment{color rgbt<1.0,0.1,0.10,0.3>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} } " <<endl;
    
    //pov<< " #declare nucleus =  texture {   pigment{color rgbt<0.50,0.1,0.10,0.81>   }   finish { diffuse 0.6    ambient 3.0 } } " <<endl;
    
    // pov<< "#declare chromatins = texture {    pigment{color rgbt<0.0,0.10,1.0,0.0>   }    finish {   diffuse 10.0  phong 2.1 phong_size 100  ambient 1.8     } }" <<endl;
    
    
    // pov<< " #declare ECM = texture { pigment{color rgbt<10/255,10/255,10/255,0>  }   }" <<endl;
    
    
    //pov<< " #declare chromatins =texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }  " <<endl;
    
    pov<< " #declare membranebonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
    // pov<< " #declare nucleusbonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
    //  pov<< "  #declare ECMbonds = texture {    pigment{color rgb<50/255,50/255,50/255>  }}  " <<endl;
    
    // pov<< "#declare  radiusChromatin=" << Chromatin_Scaling_Factor*sigmachromatin*0.5 <<";"<<endl;
    pov<< "#declare  radiBondMeM=0.03"  <<";"<<endl;
    pov<< "#declare  radiBondECM=0.1"  <<";"<<endl;
    
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    
    ///=============MEMBRANE===============
    
    pov << "union {   // MEMBRANE____________________________________________________________________" <<endl;
    
    for(int i=0;i<Num_of_Triangles;i++)
    {
        
        pov<< " triangle { < "<<Node_Position[Triangle_list[i][0]][0]<<","<<Node_Position[Triangle_list[i][0]][1] << ","<<Node_Position[Triangle_list[i][0]][2]<< ">,";
        pov<< " <"<<Node_Position[Triangle_list[i][1]][0]<<","<<Node_Position[Triangle_list[i][1]][1] << ","<<Node_Position[Triangle_list[i][1]][2]<< ">,";
        pov<< " <"<<Node_Position[Triangle_list[i][2]][0]<<","<<Node_Position[Triangle_list[i][2]][1] << ","<<Node_Position[Triangle_list[i][2]][2]<< ">";
        pov<< " } " <<endl;
    }
    pov<< " texture {membrane} no_shadow " <<endl;
    pov << " }" <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << "union { ///=========================================membrane bonds======================================="<<endl;
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        
        w1=(int)  Node_Bond_list[i][0];
        w2=(int)  Node_Bond_list[i][1];
        
        pov<< "cylinder { <" << Node_Position[w1][0]<< ","<< Node_Position[w1][1]<< "," << Node_Position[w1][2]<< "> ," ;
        pov<< "  <" << Node_Position[w2][0]<< ","<< Node_Position[w2][1]<< ","<< Node_Position[w2][2]<<"> ," <<"radiBondMeM";
        pov<< "} " <<endl;
        
        
    }
    
    pov<< " texture {membranebonds} no_shadow " <<endl;
    pov << " }" <<endl;
}



void Membrane::export_for_resume(int MD_step){
    ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Membrane_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Triangles<<endl;
    for (int i=0; i<Num_of_Triangles; i++) {
        write_resume_file<<Triangle_list[i][0]<<"\t"<<Triangle_list[i][1]<<"\t"<<Triangle_list[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    //In the import function we should call the neighbour list constructor
    write_resume_file<<Num_of_Triangle_Pairs<<endl;
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        write_resume_file<<Triangle_pair_list[i][0]<<"\t"<<Triangle_pair_list[i][1]<<"\n";
        write_resume_file<<Triangle_Pair_Nodes[i][0]<<"\t"<<Triangle_Pair_Nodes[i][1]<<"\t"<<Triangle_Pair_Nodes[i][2]<<"\t"<<Triangle_Pair_Nodes[i][3]<<"\t"<<"\n";
    }
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}

void Membrane::export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count){
    ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Membrane_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    for (int i=atom_count; i<atom_count+Num_of_Nodes; i++) {
        Node_Position[i-atom_count][0] = atoms[i].posInNm[0];
        Node_Position[i-atom_count][1] = atoms[i].posInNm[1];
        Node_Position[i-atom_count][2] = atoms[i].posInNm[2];
        
        Node_Velocity[i-atom_count][0] = atoms[i].velocityInNmperPs[0];
        Node_Velocity[i-atom_count][1] = atoms[i].velocityInNmperPs[1];
        Node_Velocity[i-atom_count][2] = atoms[i].velocityInNmperPs[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Triangles<<endl;
    for (int i=0; i<Num_of_Triangles; i++) {
        write_resume_file<<Triangle_list[i][0]<<"\t"<<Triangle_list[i][1]<<"\t"<<Triangle_list[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    //In the import function we should call the neighbour list constructor
    write_resume_file<<Num_of_Triangle_Pairs<<endl;
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        write_resume_file<<Triangle_pair_list[i][0]<<"\t"<<Triangle_pair_list[i][1]<<"\n";
        write_resume_file<<Triangle_Pair_Nodes[i][0]<<"\t"<<Triangle_Pair_Nodes[i][1]<<"\t"<<Triangle_Pair_Nodes[i][2]<<"\t"<<Triangle_Pair_Nodes[i][3]<<"\t"<<"\n";
    }
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}

void Membrane::write_geometrics(){
    calculate_volume_and_surface_area();
    calculate_surface_area_with_voronoi();
    string geo_file_name=generalParameters.trajectory_file_name+"_"+label+"_GeoProps.txt";
    std::ofstream wgeo;
    wgeo.open(geo_file_name.c_str(), std::ios::app);
    
    wgeo<<"NodeIndex NodeNormalVectorComponenets Nx Ny Nz NodeVoronoiAreaInNm2"<<endl;
    wgeo<<"timeInPs "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<"  Volume_(Nm^3) "<<volume<<"  Area_(Nm^2) "<<surface_area_voronoi<<endl;;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        wgeo<<i<<"\t"<<node_voronoi_normal_vec[i][0]<<"\t"<<node_voronoi_normal_vec[i][1]<<"\t"<<node_voronoi_normal_vec[i][2]<<"\t"<<node_voronoi_area[i]<<endl;
    }
}

void Membrane::Export_Mechanical_Energy(string filename, int frame){
    
    FILE* pFile;
    pFile = fopen (filename.c_str(),"a");
    //    fprintf(pFile,"time ps; BendingEnergy BondEnergy TotalEnergy KJ/mole\n");
    
    
    fprintf(pFile,"%8.3f\t%8.3f\t%8.3f\t%8.3f\n",
            analysis_coord_frames_time[frame],
            Total_Bending_Energy,
            Total_Potential_Energy,
            Total_Potential_Energy);
    
    fclose (pFile);
    
}


void Membrane::export_mesh_properties(string Mesh_file_name){
    string write_name = Mesh_file_name;
    write_name.erase(write_name.end()-4,write_name.end());
    write_name+="_properties.txt";
    ofstream write_prop(write_name.c_str());
    write_prop<<"# max node order "<<nodeOrder_NodeIndex_NodeNeighbourList.size()-1<<endl;
//    for (int i=0; i<nodeOrder_NodeIndex_NodeNeighbourList.size(); i++) {
//        write_prop<<i<<" "<<nodeOrder_NodeIndex_NodeNeighbourList[i].size()<<endl;
//    }
    for (int i=0; i<nodeOrder_NodeIndex_NodeNeighbourList.size(); i++) {
        if (nodeOrder_NodeIndex_NodeNeighbourList[i].size()==0) {
            continue;
        } else {
            write_prop<<"# "<<nodeOrder_NodeIndex_NodeNeighbourList[i].size()<<" nodes with order "<<to_string(i)<<"; First int: node index; and the rest are the neighbours"<<endl;
            for (int j=0; j<nodeOrder_NodeIndex_NodeNeighbourList[i].size(); j++) {
                for (int k=0; k<nodeOrder_NodeIndex_NodeNeighbourList[i][j].size(); k++) {
                    write_prop<<nodeOrder_NodeIndex_NodeNeighbourList[i][j][k]<<" ";
                }
                write_prop<<endl;
            }
        }
    }
    
    
    cout<<write_name<<endl;exit(0);
}
