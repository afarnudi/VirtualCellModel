#include <iostream>
#include <fstream>
//#include <stdlib.h>
//#include <string>
//#include <iostream>
//#include <fstream>
#include <vector>

//#include <math.h>
//#include <map>
//#include <iomanip>
//#include <iterator>
#include "General_constants.h"
#include "Membrane.h"
#include "Actin.h"
#include "ECM.h"
#include "Chromatin.h"
#include "Configfile.hpp"

using namespace std ;

string get_configfile_name(void);
bool   include_comments(void);


void configfile_generator(int status){
    
    string genfilename = get_configfile_name();
    
    
    ofstream write_configs;
    write_configs.open(genfilename.c_str());
    if (status==0) {
        bool comments      = include_comments();
        write_configs<<"#This is the VCM configuration template."<<endl;
        write_configs<<"#Whatever following \"#\" will be ignored by the programme."<<endl;
        write_configs<<"#You may use a \"word\", number, or character to label class instances. Note that this lable will only be used by the programme within the configuration file enviroment to distiguish between classes.  If no label is provided, the programme will automatically assign a number label beginning from 0."<<endl<<endl;
        
        write_configs<<"-GeneralParameters"<<endl<<endl;
        GeneralParameters defaultparams;
        for (auto & element : defaultparams.insertOrder) {
            if (comments) {
                write_configs<<defaultparams.GenParams[element][1]<<endl;
            }
            write_configs<<element<<" "<<defaultparams.GenParams[element][0]<<endl;
            write_configs<<endl;
        }
        
        
        write_configs<<"-Membrane"<<endl<<endl;
        Membrane mem;
        auto map = mem.get_map();
        auto order = mem.get_insertOrder();
        for (auto & element : order) {
            if (comments) {
                write_configs<<map[element][1]<<endl;
            }
            write_configs<<element<<" "<<map[element][0]<<endl;
            write_configs<<endl;
        }
        write_configs<<"-Actin"<<endl<<endl;
        Actin act;
        map = act.get_map();
        order = act.get_insertOrder();
        for (auto & element : order) {
            if (comments) {
                write_configs<<map[element][1]<<endl;
            }
            write_configs<<element<<" "<<map[element][0]<<endl;
            write_configs<<endl;
        }
        
        write_configs<<"-ECM"<<endl<<endl;
        ECM ecm;
        map = ecm.get_map();
        order = ecm.get_insertOrder();
        for (auto & element : order) {
            if (comments) {
                write_configs<<map[element][1]<<endl;
            }
            write_configs<<element<<" "<<map[element][0]<<endl;
            write_configs<<endl;
        }
        
        write_configs<<"-Chromatin"<<endl<<endl;
        Chromatin chromo;
        map = chromo.get_map();
        order = chromo.get_insertOrder();
        for (auto & element : order) {
            if (comments) {
                write_configs<<map[element][1]<<endl;
            }
            write_configs<<element<<" "<<map[element][0]<<endl;
            write_configs<<endl;
        }
        
        write_configs<<"-InteractionTable"<<endl<<endl;
        if (comments) {
            write_configs<<"#Define the non-bonded interaction between the Classes. The classes in the table should be written in the following order: Membranes (M0, M1, ...), Actins (A0, A1, ...), ECMs (E0, E1, ...), Chromatins (C0, C1, ...)"<<endl;
            write_configs<<"#Non-bonded interactions:"<<endl;
            write_configs<<"#0:  No interaction."<<endl;
            write_configs<<"#LJ: Lennard Jones 12-6."<<endl;
            write_configs<<"#EV: Excluded Volume interaction (Shifted Lennard Jones 12-6)."<<endl;
            write_configs<<"#Example of an interaction map with one of each calss."<<endl;
        }
        write_configs<<"    M0 A0 E0 C0"<<endl;
        write_configs<<"M0  EV         "<<endl;
        write_configs<<"A0  0  0       "<<endl;
        write_configs<<"E0  0  LJ  0    "<<endl;
        write_configs<<"C0  EV  0  LJ  0"<<endl;
        write_configs<<endl;
        if (comments) {
            write_configs<<"#Example of an interaction map with multiple calsses."<<endl;
            write_configs<<"#    M0 M1 A0 A1 E0 C0  C1  C2  C3"<<endl;
            write_configs<<"#M0  EV         "<<endl;
            write_configs<<"#M1  EV EV        "<<endl;
            write_configs<<"#A0  0  0  0     "<<endl;
            write_configs<<"#A1  0  0  0  0   "<<endl;
            write_configs<<"#E0  0  LJ 0  0  0"<<endl;
            write_configs<<"#C0  EV EV 0  0  0  EV"<<endl;
            write_configs<<"#C1  EV EV 0  0  0  EV  EV"<<endl;
            write_configs<<"#C2  EV EV 0  0  0  EV  EV  EV"<<endl;
            write_configs<<"#C3  EV EV 0  0  0  EV  EV  EV  EV"<<endl;
            write_configs<<"#Where:"<<endl;
            write_configs<<"#Excluded volume interaction is set between all Chromatins and Membrane calsses. Also, the interaction is set for each class memeber, that is the Chromatin or Membrane nodes cannot cross themselves."<<endl;
            write_configs<<"#The ECM is set to have Lennard jones interaction with Membrane 1."<<endl;
            write_configs<<"#All other classes have no non-bonded interactions."<<endl;
            write_configs<<endl<<endl;
            write_configs<<"#Note: If you have multiple instances of one class that have almost the same configurations but differ in just a couple, you can use the inheritance feature. The inheritance feature allows you to copy configurations from one class and add/overwiwte new configrations. Inheritance can only be applied to classes of the same kind:"<<endl;
            write_configs<<"#Example:"<<endl;
            write_configs<<"#-Membrane A"<<endl;
            write_configs<<"#config1 10"<<endl;
            write_configs<<"#config2 true"<<endl;
            write_configs<<"#config3 200"<<endl;
            write_configs<<"#-Membrane B :: A"<<endl;
            write_configs<<"#config2 false"<<endl;
            write_configs<<"#The programme will copy the configurations of Membrane A to Membrane B and add/overwrite the configurations that are listed in Membrane B."<<endl;
            write_configs<<"#The programme will interpret Membrane B as follows:"<<endl;
            write_configs<<"#-Membrane B"<<endl;
            write_configs<<"#config1 10"<<endl;
            write_configs<<"#config2 false"<<endl;
            write_configs<<"#config3 200"<<endl<<endl;
        }
        
        
        exit(0);
    }
    else if (status==1)
    {
        write_configs<<"#Whatever following \"#\" will be ignored by the programme."<<endl;
        write_configs<<"#You may use a \"word\", number, or character to label class instances. Note that this lable will only be used by the programme within the configuration file enviroment to distiguish between classes.  If no label is provided, the programme will automatically assign a number label beginning from 0."<<endl<<endl;
        write_configs<<"-GeneralParameters"<<endl<<endl;
        GeneralParameters defaultparams;
        cout<<TFILE;
        string input;
        cout<<TRESET;
        cout<<"VCM will generate all output files with the following format:\n"
              "ProjectName/Date+time+index/Date+time+index+file_identifier.extension\n"
              "The Date+time will be determined at the beginning of each simulation using the machine's clock. You may use a \"ProjectName\" to customise the outputfiles."<<endl;
        cout<<"Please enter a \"ProjectName\":\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
        cin>>input;
        defaultparams.GenParams["ProjectName"][0]=input;
        cout<<TRESET;
        
        cout<<"Please set the simulation time length masured in pico seconds:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
        cin>>input;
        defaultparams.GenParams["SimulationTimeInPs"][0]=input;
        cout<<TRESET;
        
        cout<<"Please set the Integration step size measured in femto seconds:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
        cin>>input;
        defaultparams.GenParams["StepSizeInFs"][0]=input;
        cout<<TRESET;
        
        cout<<"Please set how often you want to save the trajectories, etc to disk. The time intervales are measured in femto seconds:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
        cin>>input;
        defaultparams.GenParams["ReportIntervalInFs"][0]=input;
        cout<<TRESET;
        
        
        for (int i =0; i<4; i++) {
            write_configs<<defaultparams.insertOrder[i]<<" "<<defaultparams.GenParams[defaultparams.insertOrder[i]][0]<<endl;
            write_configs<<endl;
        }
        for (int i =4; i<defaultparams.insertOrder.size(); i++) {
            write_configs<<"#"<<defaultparams.insertOrder[i]<<" "<<defaultparams.GenParams[defaultparams.insertOrder[i]][0]<<endl;
        }
        write_configs<<endl;
        
        int all_class_count=0;
        int number_of_membranes=0;
        cout<<"How many membranes do you want to simulate?\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
        cin>>input;
        number_of_membranes = stoi(input);
        all_class_count+=number_of_membranes;
        cout<<TRESET<<endl;
        
        if (number_of_membranes!=0) {
            vector<Membrane> mems;
            mems.resize(number_of_membranes);
            cout<<"Configurating Membrane 0"<<endl;
            write_configs<<"-Membrane 0"<<endl<<endl;
            
            cout<<"Please set the path to the mesh file. Supported formats: Blender's ply and Gmsh 2:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
            cin>>input;
            mems[0].Params["MeshFile"][0]=input;
            cout<<TRESET;
            
            write_configs<<mems[0].insertOrder[0]<<" "<<mems[0].Params[mems[0].insertOrder[0]][0]<<endl;
            write_configs<<endl;
            
            cout<<"Please set the bond potential (harmonic) coefficient:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
            cin>>input;
            mems[0].Params["SpringCoeff"][0]=input;
            cout<<TRESET;
            
            write_configs<<mems[0].insertOrder[6]<<" "<<mems[0].Params[mems[0].insertOrder[6]][0]<<endl;
            write_configs<<endl;
            
            cout<<"Please set the bending potential (harmonic dihedral) rigidity coefficient:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
            cin>>input;
            mems[0].Params["BendingCoeff"][0]=input;
            cout<<TRESET;
            
            write_configs<<mems[0].insertOrder[8]<<" "<<mems[0].Params[mems[0].insertOrder[8]][0]<<endl;
            write_configs<<endl;
            
            for (int i =1; i<6; i++) {
                write_configs<<"#"<<mems[0].insertOrder[i]<<" "<<mems[0].Params[mems[0].insertOrder[i]][0]<<endl;
            }
            write_configs<<"#"<<mems[0].insertOrder[7]<<" "<<mems[0].Params[mems[0].insertOrder[7]][0]<<endl;
            for (int i =9; i<mems[0].insertOrder.size(); i++) {
                write_configs<<"#"<<mems[0].insertOrder[i]<<" "<<mems[0].Params[mems[0].insertOrder[i]][0]<<endl;
            }
            write_configs<<endl;
            for (int memid=1; memid<number_of_membranes; memid++) {
                cout<<"Configurating Membrane "<<memid<<endl;
                string entree;
                bool inherit = false;
                cout<<"Do you wish to inherit configurations from another Membrane (y,n)? ";
                cin>>entree;
                cout<<TRESET;
                if (entree == "y" || entree == "yes" || entree=="Yes") {
                    inherit=true;
                    cout<<"Please enter the Memebrane id:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
                    cin>>input;
                    cout<<TRESET;
                    write_configs<<"-Membrane "<<memid<<" :: "<<stoi(input)<<endl<<endl;
                    auto map = mems[memid].get_map();
                    auto order = mems[memid].get_insertOrder();
                    for (auto & element : order) {
                        write_configs<<"#"<<element<<" "<<map[element][0]<<endl;
                    }
                    write_configs<<endl;
                } else {
                
                    cout<<"Configurating Membrane "<<memid<<endl;
                    write_configs<<"-Membrane "<<memid<<endl<<endl;
                    
                    cout<<"Please set the path to the mesh file. Supported formats: Blender's ply and Gmsh 2:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
                    cin>>input;
                    mems[0].Params["MeshFile"][0]=input;
                    cout<<TRESET;
                    
                    write_configs<<mems[0].insertOrder[0]<<" "<<mems[0].Params[mems[0].insertOrder[0]][0]<<endl;
                    write_configs<<endl;
                    
                    cout<<"Please set the bond potential (harmonic) coefficient:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
                    cin>>input;
                    mems[0].Params["SpringCoeff"][0]=input;
                    cout<<TRESET;
                    
                    write_configs<<mems[0].insertOrder[6]<<" "<<mems[0].Params[mems[0].insertOrder[6]][0]<<endl;
                    write_configs<<endl;
                    
                    cout<<"Please set the bending potential (harmonic dihedral) rigidity coefficient:\n"<<TBLINK<<TBOLD<<">> "<<TRESET<<TFILE;
                    cin>>input;
                    mems[0].Params["BendingCoeff"][0]=input;
                    cout<<TRESET;
                    
                    write_configs<<mems[0].insertOrder[8]<<" "<<mems[0].Params[mems[0].insertOrder[8]][0]<<endl;
                    write_configs<<endl;
                    
                    for (int i =1; i<6; i++) {
                        write_configs<<"#"<<mems[0].insertOrder[i]<<" "<<mems[0].Params[mems[0].insertOrder[i]][0]<<endl;
                    }
                    write_configs<<"#"<<mems[0].insertOrder[7]<<" "<<mems[0].Params[mems[0].insertOrder[7]][0]<<endl;
                    for (int i =9; i<mems[0].insertOrder.size(); i++) {
                        write_configs<<"#"<<mems[0].insertOrder[i]<<" "<<mems[0].Params[mems[0].insertOrder[i]][0]<<endl;
                    }
                    write_configs<<endl;
                }
                
            }
        }
        
        write_configs<<"-InteractionTable"<<endl<<endl;
        
        
        for (int i=0; i<all_class_count; i++) {
            if (i<number_of_membranes) {
                write_configs<<"M"<<i;
            }
            for (int j=0; j<i+1; j++) {
                write_configs<<" 0";
            }
            write_configs<<endl;
        }
        
        write_configs<<endl;
        write_configs<<"#Non-bonded interactions:"<<endl;
        write_configs<<"#0:  No interaction."<<endl;
        write_configs<<"#LJ: Lennard Jones 12-6."<<endl;
        write_configs<<"#EV: Excluded Volume interaction (Shifted Lennard Jones 12-6)."<<endl;
        
        write_configs<<endl<<endl;
        write_configs<<"#Note: If you have multiple instances of one class that have almost the same configurations but differ in just a couple, you can use the inheritance feature. The inheritance feature allows you to copy configurations from one class and add/overwiwte new configrations. Inheritance can only be applied to classes of the same kind:"<<endl;
        write_configs<<"#Example:"<<endl;
        write_configs<<"#-Membrane A"<<endl;
        write_configs<<"#config1 10"<<endl;
        write_configs<<"#config2 true"<<endl;
        write_configs<<"#config3 200"<<endl;
        write_configs<<"#-Membrane B :: A"<<endl;
        write_configs<<"#config2 false"<<endl;
        write_configs<<"#The programme will copy the configurations of Membrane A to Membrane B and add/overwrite the configurations that are listed in Membrane B."<<endl;
        write_configs<<"#The programme will interpret Membrane B as follows:"<<endl;
        write_configs<<"#-Membrane B"<<endl;
        write_configs<<"#config1 10"<<endl;
        write_configs<<"#config2 false"<<endl;
        write_configs<<"#config3 200"<<endl<<endl;
    }
    
    
}


string get_configfile_name(void){
    string genfilename;
    cout<<"Please, enter the configuration file name:\n"<<TBLINK<<TBOLD<<">> "<<TRESET;
    cout<<TFILE;
    cin>>genfilename;
    cout<<TRESET;
    if (genfilename=="") {
        cout<<"Will use default: 'MyConfigfile.txt'"<<endl;
        genfilename="c";
    }
    ifstream infile(genfilename.c_str());
    bool loop=infile.is_open();
    
    while( loop ){
        string entree;
        cout<<"'"<<TBOLD<<TWARN<<genfilename<<TRESET<<"' already exists. Do you want to overwrite it (y,n)?"<<endl;
        cout<<TFILE;
        cin>>entree;
        cout<<TRESET;
        if (entree == "y" || entree == "yes" || entree=="Yes") {
            loop=false;
        } else {
            cout<<"Please, enter another name:\n"<<TBLINK<<TBOLD<<">> "<<TRESET;
            cout<<TFILE;
            cin>>genfilename;
            cout<<TRESET;
            ifstream infiletemp(genfilename.c_str());
            loop=infiletemp.is_open();
        }
    }
    return genfilename;
}

bool   include_comments(void){
    string entree;
    bool comments = false;
    cout<<"Do you want the definitions as comments (y,n)? ";
    cin>>entree;
    if (entree == "y" || entree == "yes" || entree=="Yes") {
        comments=true;
    }
    
    if (comments) {
        return true;
    }
    return false;
}
