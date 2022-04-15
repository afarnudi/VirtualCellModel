#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;


int count_mem_triangles(Triangles* triangles,
                        string     class_label
                        );

string generate_global_surface_constraint_potential_part_1(Triangles triangle, int mem_tris);
string generate_global_surface_constraint_potential_part_2(Triangles triangle);

string generate_global_volume_constraint_potential_part_1(Triangles triangle, int mem_tris);
string generate_global_volume_constraint_potential_part_2(Triangles triangle);

string generate_global_surface_constraint_potential_unified(Triangles triangle, int mem_tris);
string generate_global_volume_constraint_potential_unified(Triangles triangle, int mem_tris);

void set_surface_volume_constraint_forces(Triangles*                                 triangles,
                                          vector<OpenMM::CustomCompoundBondForce*>  &GlobalSurfaceConstraintForces,
                                          vector<OpenMM::CustomCompoundBondForce*>  &LocalSurfaceConstraintForces,
                                          vector<OpenMM::CustomCompoundBondForce*>  &GlobalVolumeConstraintForces,
                                          vector<OpenMM::CustomCompoundBondForce*>  &LocalSurfaceWCAForces,
                                          OpenMM::System                            &system
                                          ){
    set <std::string> GSCFs_classes;
    int GSCFs_index = -1;
    
    set <std::string> GVCFs_classes;
    int GVCFs_index = -1;

    
    set <std::string> LSCFs_classes;
    int LSCFs_index = -1;
    
    set <std::string> LWCAs_classes;
    int LWCAs_index = -1;
    
    set <std::string> mem_count_classes;
    
    int mem_tris=0;
    int tri_counter=0;
    
    for (int i=0; triangles[i].end_of_list_indicator != EndOfList; ++i) {
        
        auto count_class_item = mem_count_classes.find(triangles[i].class_label);
        if (count_class_item == mem_count_classes.end()) {
            mem_count_classes.insert(triangles[i].class_label);
            mem_tris = count_mem_triangles(triangles, triangles[i].class_label);
            tri_counter+=mem_tris;
        }
        
        if (triangles[i].volume_type == potentialModelIndex.Model["GlobalConstraint"]  && triangles[i].surface_type == potentialModelIndex.Model["GlobalConstraint"]) {
            auto GVCFs_item = GVCFs_classes.find(triangles[i].class_label);
            if (GVCFs_item == GVCFs_classes.end()) {

                GVCFs_classes.insert(triangles[i].class_label);
                GVCFs_index++;

                string target_volume = "target_volume_"+triangles[i].class_label;
                string target_area = "target_area_"+triangles[i].class_label;

                string surf_potential_part1 = generate_global_surface_constraint_potential_part_1(triangles[i], mem_tris);
                string volu_potential_part1 = generate_global_volume_constraint_potential_part_1(triangles[i], mem_tris);
                string potential = surf_potential_part1 + "+" + volu_potential_part1;

                GlobalVolumeConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential));
                GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_volume, triangles[i].VolumeConstraintValue);
                GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_area,   triangles[i].SurfaceConstraintValue);
                system.addForce(GlobalVolumeConstraintForces[GVCFs_index]);
                if (generalParameters.WantForce && generalParameters.force_group_count<31) {
//                    if (temp_label=="" || temp_label != triangles[i].class_label) {
                        generalParameters.force_group_count++;
//                        temp_label = triangles[i].class_label;
                        string label = triangles[i].class_label+"_VolumeArea_"+ to_string(generalParameters.force_group_count);
                        generalParameters.force_group_label.push_back(label);
//                    }
                    GlobalVolumeConstraintForces[GVCFs_index]->setForceGroup(generalParameters.force_group_count);
                }
                
                GVCFs_index++;

                string surf_potential_part2 = generate_global_surface_constraint_potential_part_2(triangles[i]);
                string volu_potential_part2 = generate_global_volume_constraint_potential_part_2(triangles[i]);
                potential = surf_potential_part2 + "+" + volu_potential_part2;

                GlobalVolumeConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(6, potential));
                GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_volume, triangles[i].VolumeConstraintValue);
                GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_area,   triangles[i].SurfaceConstraintValue);
                system.addForce(GlobalVolumeConstraintForces[GVCFs_index]);
                if (generalParameters.WantForce && generalParameters.force_group_count<31) {
//                    generalParameters.force_group_count++;
//                    string label = triangles[i].class_label+"_VolumeArea2_"+ to_string(generalParameters.force_group_count);
//                    generalParameters.force_group_label.push_back(label);
                    GlobalVolumeConstraintForces[GVCFs_index]->setForceGroup(generalParameters.force_group_count);
                }
            }

            GlobalVolumeConstraintForces[GVCFs_index-1]->addBond(triangles[i].atoms);
            if (i < tri_counter-1) {
                vector<int> atoms;
                atoms.resize(6);
                atoms[0]=triangles[i].atoms[0];
                atoms[1]=triangles[i].atoms[1];
                atoms[2]=triangles[i].atoms[2];
                int jj=i;
                for (int j=i+1; j<tri_counter; ++j) {
                    atoms[3]=triangles[j].atoms[0];
                    atoms[4]=triangles[j].atoms[1];
                    atoms[5]=triangles[j].atoms[2];
                    GlobalVolumeConstraintForces[GVCFs_index]->addBond(atoms);
                    jj++;
                }
            }

        }
        else {
            if (triangles[i].volume_type == potentialModelIndex.Model["GlobalConstraint"]) {
                
                auto GVCFs_item = GVCFs_classes.find(triangles[i].class_label);
                if (GVCFs_item == GVCFs_classes.end()) {
                    
                    GVCFs_classes.insert(triangles[i].class_label);
                    GVCFs_index++;
                    
                    string target_volume = "target_volume_"+triangles[i].class_label;
                    string potential_part1 = generate_global_volume_constraint_potential_part_1(triangles[i], mem_tris);
                    
                    GlobalVolumeConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential_part1));
                    GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_volume, triangles[i].VolumeConstraintValue);
                    system.addForce(GlobalVolumeConstraintForces[GVCFs_index]);
                    if (generalParameters.WantForce && generalParameters.force_group_count<31) {
                        generalParameters.force_group_count++;
                        GlobalVolumeConstraintForces[GVCFs_index]->setForceGroup(generalParameters.force_group_count);
                        string label = triangles[i].class_label+"_Volume_"+ to_string(generalParameters.force_group_count);
                        generalParameters.force_group_label.push_back(label);
                        
                    }
                    GVCFs_index++;
                    
                    string potential_part2 = generate_global_volume_constraint_potential_part_2(triangles[i]);
                    
                    GlobalVolumeConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(6, potential_part2));
                    GlobalVolumeConstraintForces[GVCFs_index]->addGlobalParameter(target_volume, triangles[i].VolumeConstraintValue);
                    system.addForce(GlobalVolumeConstraintForces[GVCFs_index]);
                    
                    if (generalParameters.WantForce && generalParameters.force_group_count<31) {
//                        generalParameters.force_group_count++;
                        GlobalVolumeConstraintForces[GVCFs_index]->setForceGroup(generalParameters.force_group_count);
//                        string label = triangles[i].class_label+"_Volume2_"+ to_string(generalParameters.force_group_count);
//                        generalParameters.force_group_label.push_back(label);
                        
                    }
    
                }
                
                GlobalVolumeConstraintForces[GVCFs_index-1]->addBond(triangles[i].atoms);
                if (i < tri_counter-1) {
                    vector<int> atoms;
                    atoms.resize(6);
                    atoms[0]=triangles[i].atoms[0];
                    atoms[1]=triangles[i].atoms[1];
                    atoms[2]=triangles[i].atoms[2];
                    int jj=i;
                    for (int j=i+1; j<tri_counter; ++j) {
                        atoms[3]=triangles[j].atoms[0];
                        atoms[4]=triangles[j].atoms[1];
                        atoms[5]=triangles[j].atoms[2];
                        GlobalVolumeConstraintForces[GVCFs_index]->addBond(atoms);
                        jj++;
                    }
                }
                
            }
            if (triangles[i].surface_type == potentialModelIndex.Model["GlobalConstraint"]) {

                auto GSCFs_item = GSCFs_classes.find(triangles[i].class_label);
                if (GSCFs_item == GSCFs_classes.end()) {
                    
                    GSCFs_classes.insert(triangles[i].class_label);
                    GSCFs_index++;

                    string target_surface_area = "target_area_"+triangles[i].class_label;
                    string potential_part1 = generate_global_surface_constraint_potential_part_1(triangles[i], mem_tris);
                    
                    GlobalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential_part1));
                    GlobalSurfaceConstraintForces[GSCFs_index]->addGlobalParameter(target_surface_area, triangles[i].SurfaceConstraintValue);
                    system.addForce(GlobalSurfaceConstraintForces[GSCFs_index]);
                    if (generalParameters.WantForce && generalParameters.force_group_count<31) {
                        generalParameters.force_group_count++;
                        GlobalSurfaceConstraintForces[GSCFs_index]->setForceGroup(generalParameters.force_group_count);
                        string label = triangles[i].class_label+"_Area_"+ to_string(generalParameters.force_group_count);
                        generalParameters.force_group_label.push_back(label);
                        
                    }
                    GSCFs_index++;
                    
                    string potential_part2 = generate_global_surface_constraint_potential_part_2(triangles[i]);
                    
                    GlobalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(6, potential_part2));
                    GlobalSurfaceConstraintForces[GSCFs_index]->addGlobalParameter(target_surface_area, triangles[i].SurfaceConstraintValue);
                    system.addForce(GlobalSurfaceConstraintForces[GSCFs_index]);
                    if (generalParameters.WantForce && generalParameters.force_group_count<31) {
//                        generalParameters.force_group_count++;
                        GlobalSurfaceConstraintForces[GSCFs_index]->setForceGroup(generalParameters.force_group_count);
//                        string label = triangles[i].class_label+"_Area2_"+ to_string(generalParameters.force_group_count);
//                        generalParameters.force_group_label.push_back(label);
                        
                    }
                }
                
                GlobalSurfaceConstraintForces[GSCFs_index-1]->addBond(triangles[i].atoms);
                if (i < tri_counter-1) {
                    vector<int> atoms;
                    atoms.resize(6);
                    atoms[0]=triangles[i].atoms[0];
                    atoms[1]=triangles[i].atoms[1];
                    atoms[2]=triangles[i].atoms[2];
                    int jj=i;
                    for (int j=i+1; j<tri_counter; ++j) {
                        atoms[3]=triangles[j].atoms[0];
                        atoms[4]=triangles[j].atoms[1];
                        atoms[5]=triangles[j].atoms[2];
                        GlobalSurfaceConstraintForces[GSCFs_index]->addBond(atoms);
                        jj++;
                    }
                }
            } else if (triangles[i].surface_type == potentialModelIndex.Model["LocalConstraint"]) {
                auto LSCFs_item = LSCFs_classes.find(triangles[i].class_label);
                if (LSCFs_item == LSCFs_classes.end()) {
                    
                    LSCFs_classes.insert(triangles[i].class_label);
                    LSCFs_index++;
    //                Implament local triangle area potential:
    //                \sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0/N}(a_i-\frac{A_0}{N})^2
                    string potential = "0.5*"+to_string(triangles[i].SurfaceConstraintStiffnessinKJpermolperNm2)
                    +"*("
                    +"abs("
                    +"0.5*distance(p1,p2)*distance(p1,p3)*sin(angle(p2,p1,p3))"
                    +")"
                    +"-"+to_string(triangles[i].SurfaceConstraintValue)
                    +")^2"
                    +"/"+to_string(triangles[i].SurfaceConstraintValue);
                    
                    LocalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential));
                    
                    system.addForce(LocalSurfaceConstraintForces[LSCFs_index]);
                    
                }
                LocalSurfaceConstraintForces[LSCFs_index]->addBond(triangles[i].atoms);
            }
        }
                
        if (triangles[i].area_WCA) {
            auto LWCAs_item = LWCAs_classes.find(triangles[i].class_label);
            if (LWCAs_item == LWCAs_classes.end()) {
                
                LWCAs_classes.insert(triangles[i].class_label);
                LWCAs_index++;
//                step(area_0-area)
//                string area_0 = to_string(triangles[i].surface_WCA_min_area);
                string area_0times2 = to_string(2*triangles[i].surface_WCA_min_area);
                string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
                string area = "(distance(p1,p2)*distance(p1,p3)*sin(angle(p2,p1,p3)))";
                string potential = epsilon+"*(("+area_0times2+"/"+area+")^8- ("+area_0times2+"/"+area+")^4+0.25)*step("+to_string(triangles[i].surface_WCA_min_area*pow(2,1./4.))+"-"+area+")";
                
                
                LocalSurfaceWCAForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential));
                system.addForce(LocalSurfaceWCAForces[LWCAs_index]);
            }
            LocalSurfaceWCAForces[LWCAs_index]->addBond(triangles[i].atoms);
        }
        else if (triangles[i].surface_triangle_hight_WCA) {
            auto LWCAs_item = LWCAs_classes.find(triangles[i].class_label);
            if (LWCAs_item == LWCAs_classes.end()) {
                
                LWCAs_classes.insert(triangles[i].class_label);
                LWCAs_index++;
//                step(area_0-area)
//                string area_0 = to_string(triangles[i].surface_WCA_min_area);
                string sigma_h = to_string(triangles[i].surface_triangle_hight_WCA_min_length);
                string epsilon = to_string(4*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature) ;
                string h1 = "(distance(p1,p2)*sin(angle(p3,p1,p2)))";
                string h2 = "(distance(p3,p2)*sin(angle(p1,p2,p3)))";
                string h3 = "(distance(p1,p3)*sin(angle(p2,p3,p1)))";
                string potential = epsilon+"*(("+sigma_h+"/"+h1+")^8- ("+sigma_h+"/"+h1+")^4+0.25)*step("+to_string(triangles[i].surface_triangle_hight_WCA_min_length*pow(2,1./4.))+"-"+h1+")";
                potential += "+"+epsilon+"*(("+sigma_h+"/"+h2+")^8- ("+sigma_h+"/"+h2+")^4+0.25)*step("+to_string(triangles[i].surface_triangle_hight_WCA_min_length*pow(2,1./4.))+"-"+h2+")";
                potential += "+"+epsilon+"*(("+sigma_h+"/"+h3+")^8- ("+sigma_h+"/"+h3+")^4+0.25)*step("+to_string(triangles[i].surface_triangle_hight_WCA_min_length*pow(2,1./4.))+"-"+h3+")";
                
                
                LocalSurfaceWCAForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential));
                system.addForce(LocalSurfaceWCAForces[LWCAs_index]);
            }
            LocalSurfaceWCAForces[LWCAs_index]->addBond(triangles[i].atoms);
        }
    }
}

int count_mem_triangles(Triangles* triangles,
                        string     class_label
                        ){
    int triangle_count=0;
    for (int i=0; triangles[i].end_of_list_indicator != EndOfList; ++i) {
        if (triangles[i].class_label == class_label){
            triangle_count++;
        }
    }
    return triangle_count;
}

string generate_global_surface_constraint_potential_part_1(Triangles triangle, int mem_tris){
    //                Implament Global triangle area potential:
    //                \frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^Na_i-A_0\right]^2
    //                This implamentation will make a 'string' that is too long for OpenMM to handle.
    //                Therefor we break up the potential into 2 parts:
    //                U_{liquid}^{global}=\frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^Na_i^2+A_0^2-2A_0\sum_{i=1}^Na_i+2\sum_{i\neq j}a_ia_j\right]
    //                =\frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^N(a_i-A_0)^2-(N-1)A_0^2+2\sum_{i\neq j}a_ia_j\right]
    //                =\sum_{i=1}^N\frac{1}{2}\frac{k_A}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)+\sum_{i\neq j}\frac{k_A}{A_0}a_ia_j
    
    
    
    // \sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)
    double N_1perN = (mem_tris-1)/double(mem_tris);
    string target_surface_area = "target_area_"+triangle.class_label;
    string potential_part1 = "0.5*"+to_string(triangle.SurfaceConstraintStiffnessinKJpermolperNm2)
    +"*("//begin multiplyer
    +"("
    //begin p1
//    +"abs("
    +"0.5*distance(p1,p2)*distance(p1,p3)*abs(sin(angle(p2,p1,p3)))"
//    +")"
//    +"-"+to_string(triangle.SurfaceConstraintValue)
    +"-"+target_surface_area
    +")^2"
    //end p1
    +"-"
    //begin p2
    +to_string(N_1perN)+"*"+target_surface_area+"^2"
    //end p2
    +")"//end multiplyer
//    +"/"+to_string(triangle.SurfaceConstraintValue);
    +"/"+target_surface_area;
    
    return potential_part1;
}


string generate_global_surface_constraint_potential_part_2(Triangles triangle){
    //                Implament Global triangle area potential:
    //                \frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^Na_i-A_0\right]^2
    //                This implamentation will make a 'string' that is too long for OpenMM to handle.
    //                Therefor we break up the potential into 2 parts:
    //                U_{liquid}^{global}=\frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^Na_i^2+A_0^2-2A_0\sum_{i=1}^Na_i+2\sum_{i\neq j}a_ia_j\right]
    //                =\frac{1}{2}\frac{k_A}{A_0}\left[\sum_{i=1}^N(a_i-A_0)^2-(N-1)A_0^2+2\sum_{i\neq j}a_ia_j\right]
    //                =\sum_{i=1}^N\frac{1}{2}\frac{k_A}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)+\sum_{i\neq j}\frac{k_A}{A_0}a_ia_j
    
    
    // \sum_{i\neq j}\frac{k_b}{A_0}a_ia_j
    
    string target_surface_area = "target_area_"+triangle.class_label;
    string potential_part2 = to_string(triangle.SurfaceConstraintStiffnessinKJpermolperNm2)
    +"*("//begin multiplyer
    //begin surface of first triangle
    +"0.25"
    +"*distance(p1,p2)*distance(p1,p3)"
    +"*abs(sin(angle(p2,p1,p3)))"
    +"*distance(p4,p5)*distance(p4,p6)"
    +"*abs(sin(angle(p5,p4,p6)))"
    +")"//end multiplyer
//    +"/"+to_string(triangle.SurfaceConstraintValue);
    +"/"+target_surface_area;
    
    return potential_part2;
}

string generate_global_volume_constraint_potential_part_1(Triangles triangle, int mem_tris){
    //                Implament Global triangle area potential:
    //                \frac{1}{2}k_V\left[\sum_{i=1}^Nv_i-V_0\right]^2
    //                This implamentation will make a 'string' that is too long for OpenMM to handle.
    //                Therefor we break up the potential into 2 parts:
    //                U_{liquid}^{global}=\frac{1}{2}k_V\left[\sum_{i=1}^Nv_i^2+V_0^2-2V_0\sum_{i=1}^Nv_i+2\sum_{i\neq j}v_iv_j\right]
    //                =\frac{1}{2}\frac{k_V}{V_0}\left[\sum_{i=1}^N(v_i-V_0)^2-(N-1)V_0^2+2\sum_{i\neq j}v_iv_j\right]
    //                =\sum_{i=1}^N\frac{1}{2}k_V\left((v_i-V_0)^2-\frac{(N-1)}{N}V_0^2\right)+\sum_{i\neq j}k_Vv_iv_j
    
    
    
    // \sum_{i=1}^N\frac{1}{2}k_V\left((v_i-V_0)^2-\frac{(N-1)}{N}V_0^2\right)
    double N_1perN = (mem_tris-1)/double(mem_tris);
    
    string target_volume = "target_volume_"+triangle.class_label;
    string potential_part1 = "0.5*"+to_string(triangle.VolumeConstraintStiffnessinKJpermolperNm3)
    +"*("//begin multiplyer
    +"("
    //begin p1
    // The minus sign is because the cross product of the triangle area is setup so that it always points towards the outside of the membrane.
    +"-(1./6)"
    +"*(x1*(y3*z2-z3*y2)+x2*(y1*z3-z1*y3)+x3*(y2*z1-z2*y1))"
    +"-"
//    +to_string(triangle.VolumeConstraintValue)
    +target_volume
    +")^2"
    //end p1
    +"-"
    //begin p2
    +to_string(N_1perN)+"*"+target_volume+"^2"
    //end p2
    +")"//end multiplyer
//    +"/"+to_string(triangle.VolumeConstraintValue);
    +"/"+target_volume;
    
    return potential_part1;
}


string generate_global_volume_constraint_potential_part_2(Triangles triangle){
    //                Implament Global triangle area potential:
    //                \frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^Na_i-A_0\right]^2
    //                This implamentation will make a 'string' that is too long for OpenMM to handle.
    //                Therefor we break up the potential into 2 parts:
    //                U_{liquid}^{global}=\frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^Na_i^2+A_0^2-2A_0\sum_{i=1}^Na_i+2\sum_{i\neq j}a_ia_j\right]
    //                =\frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^N(a_i-A_0)^2-(N-1)A_0^2+2\sum_{i\neq j}a_ia_j\right]
    //                =\sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)+\sum_{i\neq j}\frac{k_b}{A_0}a_ia_j
    
    
    // \sum_{i\neq j}\frac{k_b}{A_0}a_ia_j
    string target_volume = "target_volume_"+triangle.class_label;
    string potential_part2 = to_string(triangle.VolumeConstraintStiffnessinKJpermolperNm3)
    +"*(1./36)"
    +"*(x1*(y3*z2-z3*y2)+x2*(y1*z3-z1*y3)+x3*(y2*z1-z2*y1))"
    +"*(x4*(y6*z5-z6*y5)+x5*(y4*z6-z4*y6)+x6*(y5*z4-z5*y4))"
//    +"/"+to_string(triangle.VolumeConstraintValue);
    +"/"+target_volume;
    
    return potential_part2;
}

string generate_global_surface_constraint_potential_unified(Triangles triangle, int mem_tris){
    string target_surface_area = "target_area_"+triangle.class_label;
    string potential = "0.5*"+to_string(triangle.SurfaceConstraintStiffnessinKJpermolperNm2)
    +"*("
    +"("
    //begin tri 1
    +"abs("
    +"0.5*distance(p1,p2)*distance(p1,p3)*sin(angle(p2,p1,p3))"
    +")"
//    +"-"+to_string(triangle.SurfaceConstraintValue)
    +"-"+target_surface_area+"/"+to_string(mem_tris)
    //end tri 1
    +")"
    +"*("
    //begin tri 2
    +"abs("
    +"0.5*distance(p4,p5)*distance(p4,p6)*sin(angle(p5,p4,p6))"
    +")"
//    +"-"+to_string(triangle.SurfaceConstraintValue)
    +"-"+target_surface_area+"/"+to_string(mem_tris)
    //end tri 2
    +")"
    +")/"+target_surface_area;
    return potential;
}

string generate_global_volume_constraint_potential_unified(Triangles triangle, int mem_tris){
    string target_volume = "target_volume_"+triangle.class_label;
    string potential = "0.5*"+to_string(triangle.VolumeConstraintStiffnessinKJpermolperNm3)
    +"*("
    +"("
    //begin tri 1
    +"(x1*(y3*z2-z3*y2)+x2*(y1*z3-z1*y3)+x3*(y2*z1-z2*y1))/6"
    +"+"+target_volume+"/"+to_string(mem_tris)
    //end tri 1
    +")"
    +"*("
    //begin tri 2
    +"(x4*(y6*z5-z6*y5)+x5*(y4*z6-z4*y6)+x6*(y5*z4-z5*y4))/6"
    +"+"+target_volume+"/"+to_string(mem_tris)
    //end tri 2
    +")"
    +")/"+target_volume;
    return potential;
}
