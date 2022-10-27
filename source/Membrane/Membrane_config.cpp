#include <sstream>
#include "Membrane.h"
#include "Configfile.hpp"
#include "General_constants.h"
//#include "General_functions.hpp"
using namespace std;

void Membrane::import_config(vector<string> configlines){
    
    //first line is the 'key'
    configlines.erase(configlines.begin());
    
    //replace the default values with the parameters read from the Config file
    for (int i=0; i<configlines.size(); i++) {
        if(configlines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(configlines[i], "Membrane: Config reader");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = Params.find(split[0]);
                //Only replace parameters that actually exist in the Membrane parameters and ignore anythin else;
                if (it != Params.end()) {
                    configlines[i].erase(configlines[i].begin(),configlines[i].begin()+int(split[0].size()) );
                    it->second[0] = configlines[i];
                } else {
                    cout<<TWARN<<"Note: \""<<TFILE<<split[0]<<TWARN<<"\" is not a Membrane parameter."<<TRESET<<endl;
                    cout<<"If you wish to edit the configfile, exit. If not, press any key to continue."<<endl;
                    getchar();
                    cout<<TRESET;
                }
            }
        }
    }
    //Tell the class how to inturpret the strings in the config file
    assign_parameters();
    //Check if all the parameters are consistant with the physics/Class enviroment
    consistancy_check();
    //Call the initiliser from the old code
    initialise(Mesh_file_name);
}

void Membrane::consistancy_check(){
    ifstream readmesh(Mesh_file_name.c_str());
    if (!readmesh.is_open()) {
        string errorMessage = TWARN;
        errorMessage +="Read Error: Could not read/locate mesh file '"+Mesh_file_name+"'";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    if (epsilon_LJ_12_6Stat == "kbt") {
        epsilon_LJ_12_6 = generalParameters.BoltzmannKJpermolkelvin * generalParameters.temperature;
    } else {
        epsilon_LJ_12_6 = get_double_value(epsilon_LJ_12_6Stat, parser_name, "LJepsilon", "100");
    }
    
    if (Node_Bond_Nominal_Length_stat!= "Au" && Node_Bond_Nominal_Length_stat!= "Av") {
        Node_Bond_user_defined_Nominal_Length_in_Nm = get_double_value(Node_Bond_Nominal_Length_stat, parser_name, "NominalLengthInNm", "100, Au, or Av");
    }
    
    if (Triangle_pair_angle_stat!= "Au" && Triangle_pair_angle_stat!= "Av") {
        Triangle_pair_Nominal_angle_in_degrees = get_double_value(Triangle_pair_angle_stat, parser_name, "SpontaneousTriangleBendingAngleInDegrees", "100, Au, or Av");
    }
    if (Node_radius_stat!= "Au" && Node_radius_stat!= "Av") {
        double test = get_double_value(Node_radius_stat, parser_name, "NodeRadius", "100, Au, or Av");
    }
    if (SurfaceConstraintValue_stat!= "Au") {
        double test = get_double_value(SurfaceConstraintValue_stat, parser_name, "SurfaceConstraintValue", "100, Au");
    }
    if (VolumeConstraintValue_stat!= "Au") {
        double test = get_double_value(VolumeConstraintValue_stat, parser_name, "VolumeConstraintValue", "100, Au");
    }
    
    if (dihedral_bending_coefficient==0) {
        dihedral_bending_model=potentialModelIndex.Model["None"];
    }
    if (mean_curvature_coefficient==0) {
        mean_curvature_model=potentialModelIndex.Model["None"];
    }
    if (Spring_coefficient==0) {
        spring_model=potentialModelIndex.Model["None"];
    }
    if (dihedral_bending_model != potentialModelIndex.Model["None"]) {
        UseDihedralPotential =true;
    }
    
    if (mean_curvature_model != potentialModelIndex.Model["None"]) {
        UseMeanCurvature =true;
    }
    
    if (LockOnPotential == potentialModelIndex.Model["None"]) {
        LockOn_rigidity=0;
    }
    
    if (surface_constraint_model == potentialModelIndex.Model["None"] || surface_constraint_coefficient == 0) {
        surface_constraint_model = potentialModelIndex.Model["None"];
        surface_constraint_coefficient = 0;
    }
    
    if (volume_constraint_model == potentialModelIndex.Model["None"] || volume_constraint_coefficient == 0) {
        volume_constraint_model = potentialModelIndex.Model["None"];
        volume_constraint_coefficient = 0;
    }
    
    
}

void Membrane::assign_parameters(void){
    for (auto const& it : Params){
        vector<string> split = split_and_check_for_comments(it.second[0], "Membrane: "+it.first);
        if (split.size()==0) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \""+it.first+"\". Value in configuration file";
            errorMessage+=TFILE;
            errorMessage+=it.second[0];
            errorMessage+=TWARN;
            errorMessage+=". Please check the template and try again.";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        if (it.first == "MeshFile") {
            string extension = split[0];
            extension.erase(extension.begin(),extension.begin()+extension.find_last_of('.')+1);
            if (extension == "ply" || extension == "msh" || extension == "obj") {
                mesh_format=extension;
            } else{
                string errorMessage = TWARN;
                errorMessage+="I don't understand the \""+extension+"\" format. Please use ply (Blender), msh (Gmesh 2), or obj (MeshLab).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            Mesh_file_name = split[0];
        } else if (it.first == "xyzPostAnalysisPath") {
            string extension = split[0];
            extension.erase(extension.begin(),extension.begin()+extension.find_last_of('.')+1);
            if (extension != "xyz") {
                string errorMessage = TWARN;
                errorMessage+="I don't understand the \""+extension+"\" format. Please use an xyz file with xyz extension.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            xyzPostAnalysisPath = split[0];
        } else if (it.first == "NodeMass") {
            node_global_mass = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "NodeRadius") {
            Node_radius_stat = split[0];
        } else if (it.first == "SpringModel") {
            if (split[0]=="H") {
                spring_model = potentialModelIndex.Model["Harmonic"];
            } else if (split[0]=="KG") {
                spring_model = potentialModelIndex.Model["KremerGrest"];
            } else if (split[0]=="N") {
                spring_model = potentialModelIndex.Model["None"];
            } else if (split[0]=="G") {
                spring_model = potentialModelIndex.Model["Gompper"];
            } else if (split[0]=="A") {
                spring_model = potentialModelIndex.Model["Abraham1989"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Spring Model: I don't understand the \""+split[0]+"\" Model. Available models: H (Harmonic), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SurfaceConstraint") {
            if (split[0]=="L") {
                surface_constraint_model = potentialModelIndex.Model["LocalConstraint"];
            } else if (split[0]=="G") {
                surface_constraint_model = potentialModelIndex.Model["GlobalConstraint"];
            } else if (split[0]=="GWCA") {
                surface_constraint_model = potentialModelIndex.Model["GlobalConstraint"];
                surface_WCA = true;
                if (split.size()>=2) {
                    surface_WCA_min_area=get_double_value(split[1], parser_name, it.first, "GWCA 40");
                } else{
                    string errorMessage = TWARN;
                    errorMessage+="Membrane config parser: Surface Constraint Model: For the GWCA, the minimum area sould be provided after the model name: GWCA 40";
                    errorMessage+= TRESET;
                    throw std::runtime_error(errorMessage);
                }
            } else if (split[0]=="WCA") {
                surface_constraint_model = potentialModelIndex.Model["None"];
                surface_WCA = true;
                if (split.size()>=2) {
                    surface_WCA_min_area=get_double_value(split[1], parser_name, it.first, "WCA 40");
                } else{
                    string errorMessage = TWARN;
                    errorMessage+="Membrane config parser: Surface Constraint Model: For the GWCA, the minimum area sould be provided after the model name: WCA 40";
                    errorMessage+= TRESET;
                    throw std::runtime_error(errorMessage);
                }
            } else if (split[0]=="GWCAH") {
                surface_constraint_model = potentialModelIndex.Model["GlobalConstraint"];
                surface_triangle_hight_WCA = true;
                if (split.size()>=2) {
                    surface_triangle_hight_WCA_min_length=get_double_value(split[1], parser_name, it.first, "GWCAH 40");
                } else{
                    string errorMessage = TWARN;
                    errorMessage+="Membrane config parser: Surface Constraint Model: For the GWCA, the minimum area sould be provided after the model name: GWCAH 40";
                    errorMessage+= TRESET;
                    throw std::runtime_error(errorMessage);
                }
            } else if (split[0]=="WCAH") {
                surface_constraint_model = potentialModelIndex.Model["None"];
                surface_triangle_hight_WCA = true;
                if (split.size()>=2) {
                    surface_triangle_hight_WCA_min_length=get_double_value(split[1], parser_name, it.first, "WCAH 40");
                } else{
                    string errorMessage = TWARN;
                    errorMessage+="Membrane config parser: Surface Constraint Model: For the GWCA, the minimum area sould be provided after the model name: WCAH 40";
                    errorMessage+= TRESET;
                    throw std::runtime_error(errorMessage);
                }
            } else if (split[0]=="N") {
                surface_constraint_model = potentialModelIndex.Model["None"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Surface Constraint Model: I don't understand the \""+split[0]+"\" Model. Available models: G (global), L (local), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "VolumeConstraint") {
            if (split[0]=="G") {
                volume_constraint_model = potentialModelIndex.Model["GlobalConstraint"];
            } else if (split[0]=="N") {
                volume_constraint_model = potentialModelIndex.Model["None"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Volume Constraint Model: I don't understand the \""+split[0]+"\" Model. Available models: G (global), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SpringCoeff") {
            Spring_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "DampingCoeff") {
            Damping_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "DihedralBendingModel") {
            if (split[0]=="N") {
                dihedral_bending_model = potentialModelIndex.Model["None"];
            } else if (split[0]=="cosine") {
                dihedral_bending_model = potentialModelIndex.Model["Dihedral"];
            } else if (split[0]=="exp") {
                dihedral_bending_model = potentialModelIndex.Model["ExpDihedral"];
            } else if (split[0]=="smoothExp") {
                dihedral_bending_model = potentialModelIndex.Model["SmoothEXP"];
            } else if (split[0]=="smoothExp46") {
                dihedral_bending_model = potentialModelIndex.Model["SmoothEXP46"];
            } else if (split[0]=="smoothTheta4") {
                dihedral_bending_model = potentialModelIndex.Model["SmoothTheta4"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Dihedral Bending Model: I don't understand the \""+split[0]+"\" Model. Available models: cosine (cosine of dihedral angle), N (None), exp, and smoothExp.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "MeanCurvatureModel") {
            if (split[0]=="N") {
                mean_curvature_model = potentialModelIndex.Model["None"];
            } else if (split[0]=="Julicher1996") {
                mean_curvature_model = potentialModelIndex.Model["Julicher1996"];
            } else if (split[0]=="Espiru1987") {
                mean_curvature_model = potentialModelIndex.Model["Espiru1987"];
            } else if (split[0]=="Itzykson1986") {
                mean_curvature_model = potentialModelIndex.Model["Itzykson1986"];
            } else if (split[0]=="ItzyksonBarycentric") {
                mean_curvature_model = potentialModelIndex.Model["ItzyksonBarycentric"];
            } else if (split[0]=="JulicherVoronoi") {
                mean_curvature_model = potentialModelIndex.Model["JulicherVoronoi"];
            }
            else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Bending Model: I don't understand the \""+split[0]+"\" Model. Available models: cosine (cosine of dihedral angle), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "DihedralBendingCoeff") {
            dihedral_bending_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "MeanCurvatureCoeff") {
            mean_curvature_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "MeanspontaneousCurvature") {
            mean_spontaneous_curvature = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "SurfaceConstraintCoeff") {
            surface_constraint_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "VolumeConstraintCoeff") {
            volume_constraint_coefficient = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "NominalLengthInNm") {
            Node_Bond_Nominal_Length_stat = split[0];
        } else if (it.first == "SpontaneousTriangleBendingAngleInDegrees") {
            Triangle_pair_angle_stat = split[0];
        } else if (it.first == "SurfaceConstraintValue") {
            SurfaceConstraintValue_stat = split[0];
        } else if (it.first == "VolumeConstraintValue") {
            VolumeConstraintValue_stat = split[0];
        } else if (it.first == "VolumeConstraintRatio") {
            VolumeConstraintRatio = get_double_value(split[0], parser_name, it.first, "1");
        } else if (it.first == "ReducedNominalLengthRatio") {
            ReducedNominalLengthRatio = get_double_value(split[0], parser_name, it.first, "1");
        } else if (it.first == "SurfaceConstraintRatio") {
            SurfaceConstraintRatio = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "LinearReducedSrfaceVolume") {
            LinearReducedSrfaceVolume = get_int_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "LinearReducedNominalLength") {
            if (split.size()!=2) {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: LinearReducedNominalLength: Please specify two numbers for the initial and the final time step (in Pico seconds). ";
                errorMessage+=TFILE+to_string(split.size())+TWARN+" argument(s) provided:\n";
                for (int split_id=0; split_id<split.size(); split_id++) {
                    errorMessage+=TFILE+split[split_id]+TWARN+", ";
                }
                errorMessage+="\n";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            LinearReducedNominalLength[0] = get_int_value(split[0], parser_name, it.first, "100");
            LinearReducedNominalLength[1] = get_int_value(split[1], parser_name, it.first, "100");
        } else if (it.first == "xyzPostAnalysisPathFrame") {
            if (split[0]=="L") {
                xyzPostAnalysisPathFrame=-1;
            } else {
                xyzPostAnalysisPathFrame = get_int_value(split[0], parser_name, it.first, "2");
            }
        } else if (it.first == "XYZinMembrane") {
            UseXYZinMembrane = true;
            XYZinMembrane.resize(3);
            XYZinMembrane[0] = get_double_value(split[0], parser_name, it.first, "100");
            XYZinMembrane[1] = get_double_value(split[1], parser_name, it.first, "100");
            XYZinMembrane[2] = get_double_value(split[2], parser_name, it.first, "100");
        } else if (it.first == "XYZscale") {
            X_scale = get_double_value(split[0], parser_name, it.first, "100");
            Y_scale = get_double_value(split[1], parser_name, it.first, "100");
            Z_scale = get_double_value(split[2], parser_name, it.first, "100");
        } else if (it.first == "Scale") {
            rescale_factor = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "LJepsilon") {
            epsilon_LJ_12_6Stat  = split[0];
        } else if (it.first == "UpdateRadius") {
            New_Radius = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "UpdateBeginTimeInPs") {
            Begin_update_time_in_Ps = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "UpdateEndTimeInPs") {
            End_update_time_in_Ps = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "InitRandomRotation") {
            initial_random_rotation_coordinates=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "CalculateBending") {
            CalculateBending=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "Sticky") {
            sticky=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "WantGeometricProps") {
            WantGeometricProps=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "InflateMembrane") {
            InflateMembrane=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "FreezeSubLattice") {
            freezeSubLattice=get_bool_value(split[0], parser_name, it.first);
        } else if (it.first == "ExtForceModel") {
            ext_force_model = get_int_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "LockOnRigidity") {
            LockOn_rigidity = get_double_value(split[0], parser_name, it.first, "100");
        } else if (it.first == "LockOnPotential") {
            if (split[0]=="N") {
                LockOnPotential = potentialModelIndex.Model["None"];
            } else if (split[0]=="Sphere") {
                LockOnPotential = potentialModelIndex.Model["LockOnSphere"];
            } else if (split[0]=="Ellipsoid") {
                LockOnPotential = potentialModelIndex.Model["LockOnEllipsoid"];
            } else if (split[0]=="ULM2_0") {
                LockOnPotential = potentialModelIndex.Model["LockOnULM2_0"];
                if (split.size()>1) {
                    LockOnULMU = get_double_value(split[1], parser_name, it.first, "1");
                } else {
                    LockOnULMU = 1;
                }
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: LockOnPotential: I don't understand the \""+split[0]+"\" potential. Available potentials: Sphere (for a spherical surface potential), Ellipsoid, N (None), and ULM2_0 followed by the value of u (default 1).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "LockOnCoordinate") {
            LockOnCoordinate.resize(3,0);
            LockOnCoordinate[0] = get_double_value(split[0], parser_name, it.first, "100");
            LockOnCoordinate[1] = get_double_value(split[1], parser_name, it.first, "100");
            LockOnCoordinate[2] = get_double_value(split[2], parser_name, it.first, "100");
            
        } else if (it.first == "ExtForceRigidity") {
            ext_force_rigidity.resize(3,0);
            ext_force_rigidity[0] = get_double_value(split[0], parser_name, it.first, "100");
            ext_force_rigidity[1] = get_double_value(split[1], parser_name, it.first, "100");
            ext_force_rigidity[2] = get_double_value(split[2], parser_name, it.first, "100");
            
        } else if (it.first == "VelocityShiftVector") {
            Shift_velocities_xyzVector.resize(3,0);
            Shift_velocities_xyzVector[0] = get_double_value(split[0], parser_name, it.first, "100");
            Shift_velocities_xyzVector[1] = get_double_value(split[1], parser_name, it.first, "100");
            Shift_velocities_xyzVector[2] = get_double_value(split[2], parser_name, it.first, "100");
        } else if (it.first == "CoordinateTranslateVector") {
            Shift_position_xyzVector.resize(3,0);
            Shift_position_xyzVector[0] = get_double_value(split[0], parser_name, it.first, "100");
            Shift_position_xyzVector[1] = get_double_value(split[1], parser_name, it.first, "100");
            Shift_position_xyzVector[2] = get_double_value(split[2], parser_name, it.first, "100");
        } else if (it.first == "NUL") {
            if (get_int_value(split[0], parser_name, it.first, "100")!=0) {
                AddRandomModes=true;
                NumberOfRandomModes = get_int_value(split[0], parser_name, it.first, "100");
                UlmOfRandomModes    = get_double_value(split[1], parser_name, it.first, "100");
                EllMaxOfRandomModes = get_int_value(split[2], parser_name, it.first, "100");
            }
            
        } else if (it.first == "ULM") {
            if (get_int_value(split[1], parser_name, it.first, "2")!=0) {
                AddSphericalHarmonicsMode=true;
                AmplitudeOfGeneratedMode = get_double_value(split[0], parser_name, it.first, "1");
                EllOfGeneratedMode = get_int_value(split[1], parser_name, it.first, "2");
                MOfGeneratedMode = get_int_value(split[2], parser_name, it.first, "2");
            }
            
        }  else if (it.first == "CalculateGeometryBeforeModeDisconfiguration") {
            CalculateGeometryBeforeModeDisconfiguration=get_bool_value(split[0], parser_name, it.first);
        }
        
    }
}
