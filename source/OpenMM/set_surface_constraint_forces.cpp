#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

int count_mem_nodes(Triangles* triangles,
                    string     class_label
                    );

int count_mem_triangles(Triangles* triangles,
                        string     class_label
                        );

string generate_surface_constraint_potential(Triangles* triangles,
                                             string     class_label,
                                             int        triangle_count
                                             );
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
    set <std::string> LSCFs_classes;
    int LSCFs_index = -1;
    
    int mem_nodes=0;
    int atom_counter=0;
    for (int i=0; triangles[i].type != EndOfList; ++i) {
        
        
        if (triangles[i].type == potentialModelIndex.Model["GlobalConstraint"]) {
            auto GSCFs_item = GSCFs_classes.find(triangles[i].class_label);
            if (GSCFs_item == GSCFs_classes.end()) {
                
                GSCFs_classes.insert(triangles[i].class_label);
                GSCFs_index++;
                
                string potential = generate_surface_constraint_potential(triangles, triangles[i].class_label, i);
//                cout<<potential<<endl;exit(0);
                mem_nodes = count_mem_nodes(triangles, triangles[i].class_label);
//                cout<<"mem_nodes = "<<mem_nodes<<endl;exit(0);
                
                GlobalSurfaceConstraintForces.push_back(new OpenMM::CustomCompoundBondForce(mem_nodes, potential));
                
                system.addForce(GlobalSurfaceConstraintForces[GSCFs_index]);
                vector<int> particles = get_global_surface_constraint_particles(mem_nodes, atom_counter);
                GlobalSurfaceConstraintForces[GSCFs_index]->addBond(particles);
                atom_counter+=mem_nodes;
            }
            
        } else if (triangles[i].type == potentialModelIndex.Model["LocalConstraint"]) {
            auto LSCFs_item = LSCFs_classes.find(triangles[i].class_label);
            if (LSCFs_item == LSCFs_classes.end()) {
                
                LSCFs_classes.insert(triangles[i].class_label);
                LSCFs_index++;
                
                string potential = "0.5*"+to_string(triangles[i].ConstraintStiffnessinKJpermolperNm2)
                    +"*("
                    +"abs("
                    +"distance(p1,p2)*distance(p1,p3)*sin(angle(p2,p1,p3))"
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


string generate_surface_constraint_potential(Triangles* triangles,
                                             string     class_label,
                                             int        triangle_count
                                             ){
    
    string area="0.5*(";
    int mem_nodes = count_mem_nodes(triangles, class_label);
    int mem_tris  = count_mem_triangles(triangles, class_label);
    for (int i=0; i<mem_tris; ++i) {
        area+="abs(";
        area+="distance(p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[1]+1) + ")*";
        area+="distance(p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[2]+1) + ")*" ;
        area+="sin(angle(p" + to_string(triangles[i].atoms[1]+1) + ",p" + to_string(triangles[i].atoms[0]+1) + ",p" + to_string(triangles[i].atoms[2]+1) + "))";
        area+=")+";
    }
    //extra '+' at the end
    area.pop_back();
    area+=")";
    string potential="0.5*";
    potential+=to_string(triangles[triangle_count].ConstraintStiffnessinKJpermolperNm2)+"*(";
    potential+=area+"-"+to_string(triangles[triangle_count].SurfaceConstraintValue)+")^2";
    potential+="/"+to_string(triangles[triangle_count].SurfaceConstraintValue);
//    cout<<potential<<endl;
    
    
//    string area="";
////    area+="0.5*(";
//    int mem_nodes = count_mem_nodes(triangles, class_label);
//    int mem_tris  = count_mem_triangles(triangles, class_label);
//    string defs="";
////    for (int i=0; i<mem_tris; ++i) {
//    for (int i=0; i<3; ++i) {
//        string p1=to_string(triangles[i].atoms[0]+1);
//        string p2=to_string(triangles[i].atoms[1]+1);
//        string p3=to_string(triangles[i].atoms[2]+1);
//        
//        string r12 = "r" + p1 + "_" + p2;
//        string r13 = "r" + p1 + "_" + p3;
//        string t213= "t" + p2 + "_" + p1 + "_" + p3;
//        
//        string r12def = r12 +"=distance(p" + p1 + ",p" + p2 + ");";
//        string r13def = r13 +"=distance(p" + p1 + ",p" + p3 + ");";
//        string t213def= t213 +"=angle(p" + p2 + ",p" + p1 + ",p" + p3 + ");";
//        defs+= r12def + r13def + t213def;
//        
////        area+="abs(";
//        area+= r12 + "*";
//        area+= r13 + "*" ;
//        area+="sin(" + t213 + ")";
////        area+=")";
//        area+="+";
//    }
//    //extra '+' at the end
//    area.pop_back();
////    area+=" )";
//    string potential="";
//    potential+="0.25*";
//    potential+=to_string(triangles[triangle_count].ConstraintStiffnessinKJpermolperNm2)+"*(";
//    potential+=area+"-"+to_string(triangles[triangle_count].SurfaceConstraintValue)+")^2";
//    potential+="/"+to_string(triangles[triangle_count].SurfaceConstraintValue);
//    potential+=";"+defs;
////    cout<<potential<<endl;
    return potential;
}