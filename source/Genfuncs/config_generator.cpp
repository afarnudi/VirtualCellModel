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
    bool comments      = include_comments();
    
    ofstream write_configs;
    write_configs.open(genfilename.c_str());
    if (status==0) {
        write_configs<<"#This is the VCM configuration template."<<endl;
        write_configs<<"#Whatever following \"\#\" will be ignored by the programme."<<endl;
        write_configs<<"#You may use a \"word\", number, or character to label class instances. Note that this lable will only be used by the programme within the configuration file enviroment to distiguish between classes.  If no label is provided, the programme will automatically assign a number label beginning from 0."<<endl;
        write_configs<<"#Note: If you have multiple instances of one class with only a couple of deviations, you can inherit the configrations from another class of the same kind:"<<endl;
        write_configs<<"#Example:"<<endl;
        write_configs<<"#-Membrane A"<<endl;
        write_configs<<"#list of configurations"<<endl;
        write_configs<<"#-Membrane B :: A"<<endl;
        write_configs<<"#Settings that differ"<<endl;
        write_configs<<"#The programme will copy the configurations of Membrane A to Membrane B and add/overwrite the configurations that are listed in Membrane B."<<endl<<endl;
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
        }
        
        
        exit(0);
    }
    else if (status==1)
    {
        
        
        
    }
    
    else if (status==2)
    {
        
    }
    
    else if (status==3)
    {
        
    }
    
}


string get_configfile_name(void){
    string genfilename;
    cout<<"Please, enter the configuration file name:\n"<<TBLINK<<TBOLD<<">> "<<TRESET;
    cin>>genfilename;
    if (genfilename=="") {
        cout<<"Will use default: 'MyConfigfile.txt'"<<endl;
        genfilename="c";
    }
    ifstream infile(genfilename.c_str());
    bool loop=infile.is_open();
    
    while( loop ){
        string entree;
        cout<<"'"<<TBOLD<<TWARN<<genfilename<<TRESET<<"' already exists. Do you want to overwrite it (y,n)?"<<endl;
        cin>>entree;
        if (entree == "y" || entree == "yes" || entree=="Yes") {
            loop=false;
        } else {
            cout<<"Please, enter another name:\n"<<TBLINK<<TBOLD<<">> "<<TRESET;
            cin>>genfilename;
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
