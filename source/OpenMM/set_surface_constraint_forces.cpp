#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

int count_mem_nodes(Triangles* triangles,
                    string     class_label
                    );

int count_mem_triangles(Triangles* triangles,
                        string     class_label
                        );

//string generate_surface_constraint_potential_part1(Triangles* triangles,
//                                                   string     class_label,
//                                                   int        triangle_count
//                                                   );
vector<int> get_global_surface_constraint_particles(int mem_nodes,
                                                    int atom_count);

void set_surface_constraint_forces(Triangles*                                 triangles,
                                   vector<OpenMM::CustomCompoundBondForce*>  &GlobalSurfaceConstraintForces,
                                   vector<OpenMM::CustomCompoundBondForce*>  &LocalSurfaceConstraintForces,
                                   OpenMM::System                            &system
                                   ){
    //DFs = DihedralForces
    set <std::string> GSCFs_classes;
    int GSCFs_index = -1;
//    set <std::string> GSCFsp2_classes;
//    int GSCFsp2_index = -1;
    
    set <std::string> LSCFs_classes;
    int LSCFs_index = -1;
    
    set <std::string> mem_count_classes;
    
    int mem_nodes=0;
    int mem_tris=0;
    int tri_counter=0;
    
    for (int i=0; triangles[i].type != EndOfList; ++i) {
        
        auto count_class_item = mem_count_classes.find(triangles[i].class_label);
        if (count_class_item == mem_count_classes.end()) {
            mem_count_classes.insert(triangles[i].class_label);
            mem_tris = count_mem_triangles(triangles, triangles[i].class_label);
            tri_counter+=mem_tris;
        }
        if (triangles[i].type == potentialModelIndex.Model["GlobalConstraint"]) {
            auto GSCFs_item = GSCFs_classes.find(triangles[i].class_label);
            if (GSCFs_item == GSCFs_classes.end()) {
                
                GSCFs_classes.insert(triangles[i].class_label);
                GSCFs_index++;
//                Implament Global triangle area potential:
//                \frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^Na_i-A_0\right]^2
//                This implamentation will make a 'string' that is too long for OpenMM to handle.
//                Therefor we break up the potential into 2 parts:
//                U_{liquid}^{global}=\frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^Na_i^2+A_0^2-2A_0\sum_{i=1}^Na_i+2\sum_{i\neq j}a_ia_j\right]
//                =\frac{1}{2}\frac{k_b}{A_0}\left[\sum_{i=1}^N(a_i-A_0)^2-(N-1)A_0^2+2\sum_{i\neq j}a_ia_j\right]
//                =\sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)+\sum_{i\neq j}\frac{k_b}{A_0}a_ia_j
                
                
//                \sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0}\left((a_i-A_0)^2-\frac{(N-1)}{N}A_0^2\right)
                double N_1timesA2 = (mem_tris-1)*triangles[i].SurfaceConstraintValue*triangles[i].SurfaceConstraintValue/double(mem_tris);
                string potential_part1 = "0.5*"+to_string(triangles[i].ConstraintStiffnessinKJpermolperNm2)
                +"*("//begin multiplyer
                +"("
                //begin p1
                +"abs("
                +"0.5*distance(p1,p2)*distance(p1,p3)*sin(angle(p2,p1,p3))"
                +")"
                +"-"+to_string(triangles[i].SurfaceConstraintValue)
                +")^2"
                //end p1
                +"-"
                //begin p2
                +to_string(N_1timesA2)
                //end p2
                +")"//end multiplyer
                +"/"+to_string(triangles[i].SurfaceConstraintValue);
                
                GlobalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(3, potential_part1));
                system.addForce(GlobalSurfaceConstraintForces[GSCFs_index]);
//                cout<<potential_part1<<endl<<endl;
                GSCFs_index++;
                
//                \sum_{i\neq j}\frac{k_b}{A_0}a_ia_j
                string potential_part2 = to_string(triangles[i].ConstraintStiffnessinKJpermolperNm2)
                +"*("//begin multiplyer
                //begin surface of first triangle
                +"0.25"
                +"*distance(p1,p2)*distance(p1,p3)"
                +"*abs(sin(angle(p2,p1,p3)))"
                +"*distance(p4,p5)*distance(p4,p6)"
                +"*abs(sin(angle(p5,p4,p6)))"
                +")"//end multiplyer
                +"/"+to_string(triangles[i].SurfaceConstraintValue);
                
                GlobalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(6, potential_part2));
                system.addForce(GlobalSurfaceConstraintForces[GSCFs_index]);
//                cout<<potential_part2<<endl;exit(0);
            }
            GlobalSurfaceConstraintForces[GSCFs_index-1]->addBond(triangles[i].atoms);
//            cout<<"G :"<<triangles[i].atoms[0]<<" "<<triangles[i].atoms[1]<<" "<<triangles[i].atoms[2]<<"\n";
            if (i == tri_counter-1) {
                continue;
            } else {
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
//                for (int k=0; k<6; k++) {
//                    cout<<atoms[k]<<" ";
//                }
//                cout<<endl;
                
            }
        } else if (triangles[i].type == potentialModelIndex.Model["LocalConstraint"]) {
            auto LSCFs_item = LSCFs_classes.find(triangles[i].class_label);
            if (LSCFs_item == LSCFs_classes.end()) {
                
                LSCFs_classes.insert(triangles[i].class_label);
                LSCFs_index++;
//                Implament local triangle area potential:
//                \sum_{i=1}^N\frac{1}{2}\frac{k_b}{A_0/N}(a_i-\frac{A_0}{N})^2
                string potential = "0.5*"+to_string(triangles[i].ConstraintStiffnessinKJpermolperNm2)
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
//            cout<<"L :"<<triangles[i].atoms[0]<<" "<<triangles[i].atoms[1]<<" "<<triangles[i].atoms[2]<<"\n";
        }
        
    }
//    exit(0);
}

vector<int> get_global_surface_constraint_particles(int mem_nodes,
                                                    int atom_count){
    vector<int> particles;
    particles.resize(mem_nodes);
    for (int i=atom_count; i<atom_count+mem_nodes; i++) {
        particles[i-atom_count]=i;
    }
    return particles;
}

int count_mem_nodes(Triangles* triangles,
                    string     class_label
                    ){
    std::set<int> mem_set_of_nodes;
    for (int i=0; triangles[i].type != EndOfList; ++i) {
        if (triangles[i].class_label == class_label){
            for (auto &node : triangles[i].atoms){
                mem_set_of_nodes.insert(node);
            }
        }
    }
    int node_count = mem_set_of_nodes.size();
    return node_count;
}

int count_mem_triangles(Triangles* triangles,
                        string     class_label
                        ){
    int triangle_count=0;
    for (int i=0; triangles[i].type != EndOfList; ++i) {
        if (triangles[i].class_label == class_label){
            triangle_count++;
        }
    }
    return triangle_count;
}


//string generate_surface_constraint_potential_part1(Triangles* triangles,
//                                                   string     class_label,
//                                                   int        triangle_count
//                                                   ){
//
//    string area="0.5*(";
//    area+="abs(";
//    area+="distance(p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[1]+1) + ")*";
//    area+="distance(p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[2]+1) + ")*" ;
//    area+="sin(angle(p" + to_string(triangles[i].atoms[1]+1) + ",p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[2]+1) + "))";
//    area+=")+";
//
//    //extra '+' at the end
//    area.pop_back();
//    area+=")";
//    string potential="0.5*";
//    potential+=to_string(triangles[triangle_count].ConstraintStiffnessinKJpermolperNm2)+"*(";
//    potential+=area+"-"+to_string(triangles[triangle_count].SurfaceConstraintValue)+")^2";
//    potential+="/"+to_string(triangles[triangle_count].SurfaceConstraintValue);
//    return potential;
//}
