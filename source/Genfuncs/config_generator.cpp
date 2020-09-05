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
#include "maps.hpp"

using namespace std ;

string get_configfile_name(void);
bool   include_comments(void);


void configfile_generator(int status){
    
    string genfilename = get_configfile_name();
    bool comments      = include_comments();
    
    ofstream write_configs;
    write_configs.open(genfilename.c_str());
    if (status==0) {
        
        write_configs<<"-General Parameters"<<endl<<endl;
        GeneralParameters defaultparams;
        for (auto const& it : defaultparams.GenParams)
        {
            write_configs<<"#"<<it.first<<" "<<it.second[0];
            if (comments) {
                write_configs<<" "<<it.second[1];
            }
            write_configs<<endl;
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