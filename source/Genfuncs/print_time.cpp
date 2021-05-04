//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include "General_functions.hpp"
#include <boost/filesystem.hpp>

void print_time(std::string filepath, std::string hardwareReportHeader, bool printToScreen,
                clock_t tStart,
                std::chrono::time_point<std::chrono::steady_clock> chronoSteadyClockStart,
                std::chrono::time_point<std::chrono::system_clock> chronoSystemClockStart
                ){
    
    string backupfilepath = filepath +"Backup";
    boost::filesystem::copy_file(filepath, backupfilepath, boost::filesystem::copy_option::overwrite_if_exists);
    
    ofstream file(filepath.c_str());

    double sim_duration_per_sec = (double)((clock() - tStart)/CLOCKS_PER_SEC);
    std::chrono::time_point<std::chrono::steady_clock> chronoSteadyClockEnd = chrono::steady_clock::now();
    std::chrono::time_point<std::chrono::system_clock> chronoSystemClockEnd = chrono::system_clock::now();
    
    
    double sec_per_day     =60*60*24;
    double sec_per_hour    =60*60;
    double sec_per_min     =60;
    
    int days =  sim_duration_per_sec/sec_per_day;
    sim_duration_per_sec -= days * sec_per_day;
    
    int hours = sim_duration_per_sec/sec_per_hour;
    sim_duration_per_sec -= hours * sec_per_hour;
    
    int mins = sim_duration_per_sec/sec_per_min ;
    sim_duration_per_sec -= mins * sec_per_min;
    
    file<<hardwareReportHeader;
    
    if (printToScreen) {
        cout << "Wall clock time of the simulation:\n";
        if (days!=0) {
            cout << days << "\tDays\n";
        }
        if (hours!=0) {
            cout << hours << "\tHours\n";
        }
        if (mins!=0) {
            cout << mins  << "\tMinutes\n";
        }
        cout << int(sim_duration_per_sec) << "\tSeconds\n";
    }
    
    file << "Wall clock time of the simulation:\n";
    if (days!=0) {
        file << days << "\tDays\n";
    }
    if (hours!=0) {
        file << hours << "\tHours\n";
    }
    if (mins!=0) {
        file << mins  << "\tMinutes\n";
    }
    file << int(sim_duration_per_sec) << "\tSeconds\n\n";
    
    
    auto chromoSteadyClockDiff = chronoSteadyClockEnd - chronoSteadyClockStart;
    days  = int(std::chrono::duration_cast<std::chrono::hours>(chromoSteadyClockDiff).count()/24);
    hours = std::chrono::duration_cast<std::chrono::hours>(chromoSteadyClockDiff).count() - days*24;
    mins  = std::chrono::duration_cast<std::chrono::minutes>(chromoSteadyClockDiff).count();
    int secs  = std::chrono::duration_cast<std::chrono::seconds>(chromoSteadyClockDiff).count();
    
    file << "\nchrono::steady_clock runtime: \n";
    if (days!=0) {
        file << days << "\tDays\n";
    }
    if (hours!=0) {
        file << hours << "\tHours\n";
    }
    if (mins!=0) {
        file << mins - hours * 60 << "\tMinutes\n";
    }
    file << secs - mins * 60 << "\tSeconds\n\n";
    if (printToScreen) {
        cout << "\nchrono::steady_clock runtime: \n";
        if (days!=0) {
            cout << days << "\tDays\n";
        }
        if (hours!=0) {
            cout << hours << "\tHours\n";
        }
        if (mins!=0) {
            cout << mins - hours * 60 << "\tMinutes\n";
        }
        cout << secs - mins * 60 << "\tSeconds\n";
    }
    
    
    auto chromoSystemClockDiff = chronoSystemClockEnd - chronoSystemClockStart;
    days =  int(std::chrono::duration_cast<std::chrono::hours>(chromoSystemClockDiff).count()/24);
    hours = std::chrono::duration_cast<std::chrono::hours>(chromoSystemClockDiff).count()-days*24;
    mins = std::chrono::duration_cast<std::chrono::minutes>(chromoSystemClockDiff).count();
    secs = std::chrono::duration_cast<std::chrono::seconds>(chromoSystemClockDiff).count();
    
    file << "\nchrono::system_clock runtime: \n";
    if (days!=0) {
        file << days << "\tDays\n";
    }
    if (hours!=0) {
        file << hours << "\tHours\n";
    }
    if (mins!=0) {
        file << mins - hours * 60 << "\tMinutes\n";
    }
    file << secs - mins * 60 << "\tSeconds\n\n";
    if (printToScreen) {
        cout << "\nchrono::system_clock runtime: \n";
        if (days!=0) {
            cout << days << "\tHours\n";
        }
        if (hours!=0) {
            cout << hours << "\tHours\n";
        }
        if (mins!=0) {
            cout << mins - hours * 60 << "\tMinutes\n";
        }
        cout << secs - mins * 60 << "\tSeconds\n";
    }
    file.close();
}

void print_wall_clock_time(double sim_duration_per_sec, bool printToScreen){
    std::string traj_name  = generalParameters.trajectory_file_name+"_hardware_runtime.txt";
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    
    
    double sec_per_day     =60*60*24;
    double sec_per_hour    =60*60;
    double sec_per_min     =60;
    
    
    int days =  sim_duration_per_sec/sec_per_day;
    sim_duration_per_sec -= days * sec_per_day;
    
    int hours = sim_duration_per_sec/sec_per_hour;
    sim_duration_per_sec -= hours * sec_per_hour;
    
    int mins = sim_duration_per_sec/sec_per_min ;
    sim_duration_per_sec -= mins * sec_per_min;
    if (printToScreen) {
        printf("Wall clock time of the simulation:\n");
        printf("%4i\tHours,\n%4i\tMinutes,\n%4i\tSeconds\n", hours,mins,int(sim_duration_per_sec) );
    }
    
    
    fprintf(pFile,"Wall clock time of the simulation:\n");
    fprintf(pFile,"%4i\tHours,\n%4i\tMinutes,\n%4i\tSeconds\n", hours,mins,int(sim_duration_per_sec) );
    fclose (pFile);
}

//using std::chrono;
using std::cout;

void print_real_time(std::chrono::time_point<std::chrono::steady_clock> chrono_clock_start,
                     std::chrono::time_point<std::chrono::steady_clock> chrono_clock_end,
                     bool printToScreen){
    std::string hardname= generalParameters.trajectory_file_name+"_hardware_runtime.txt";
    std::ofstream write_report;
    write_report.open(hardname.c_str(),std::ios_base::app);
    
    auto chromo_clock_diff = chrono_clock_end - chrono_clock_start;
    
    int hours = std::chrono::duration_cast<std::chrono::hours>(chromo_clock_diff).count();
    int mins  = std::chrono::duration_cast<std::chrono::minutes>(chromo_clock_diff).count();
    int secs  = std::chrono::duration_cast<std::chrono::seconds>(chromo_clock_diff).count();
    write_report << "\nchrono::steady_clock runtime: \n";
    write_report << hours << "\tHours\n";
    write_report << mins - hours * 60 << "\tMinutes\n";
    write_report << secs - mins * 60 << "\tSeconds\n\n";
    if (printToScreen) {
        cout << "\nchrono::steady_clock runtime: \n";
        cout << hours << "\tHours\n";
        cout << mins - hours * 60 << "\tMinutes\n";
        cout << secs - mins * 60 << "\tSeconds\n";
    }
    
    
    
}

void print_system_time(std::chrono::time_point<std::chrono::system_clock> chrono_clock_start,
                       std::chrono::time_point<std::chrono::system_clock> chrono_clock_end,
                       bool printToScreen){
    std::string hardname= generalParameters.trajectory_file_name+"_hardware_runtime.txt";
    std::ofstream write_report;
    write_report.open(hardname.c_str(),std::ios_base::app);
    
    auto chromo_clock_diff = chrono_clock_end - chrono_clock_start;
    int secs = std::chrono::duration_cast<std::chrono::seconds>(chromo_clock_diff).count();
    int hours = std::chrono::duration_cast<std::chrono::hours>(chromo_clock_diff).count();
    int mins = std::chrono::duration_cast<std::chrono::minutes>(chromo_clock_diff).count();
    
    write_report << "\nchrono::system_clock runtime: \n";
    write_report << hours << "\tHours\n";
    write_report << mins - hours * 60 << "\tMinutes\n";
    write_report << secs - mins * 60 << "\tSeconds\n\n";
    if (printToScreen) {
        cout << "\nchrono::system_clock runtime: \n";
        cout << hours << "\tHours\n";
        cout << mins - hours * 60 << "\tMinutes\n";
        cout << secs - mins * 60 << "\tSeconds\n";
    }
}
