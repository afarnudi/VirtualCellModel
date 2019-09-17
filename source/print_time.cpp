//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include "General_functions.hpp"

void print_wall_clock_time(double sim_duration_per_sec){
    
    double sec_per_day     =60*60*24;
    double sec_per_hour    =60*60;
    double sec_per_min     =60;
    
    
    int days =  sim_duration_per_sec/sec_per_day;
    sim_duration_per_sec -= days * sec_per_day;
    
    int hours = sim_duration_per_sec/sec_per_hour;
    sim_duration_per_sec -= hours * sec_per_hour;
    
    int mins = sim_duration_per_sec/sec_per_min ;
    sim_duration_per_sec -= mins * sec_per_min;
    
    printf("Wall clock time of the simulation:\n");
    //            printf("%.2f Days,\n%.2f Hours,\n%.2f Minutes,\n%.2f Seconds\n", days,hours,mins,sim_duration_per_sec);
    printf("%4i\tHours,\n%4i\tMinutes,\n%4i\tSeconds\n", hours,mins,int(sim_duration_per_sec) );
}

//using std::chrono;
using std::cout;

void print_real_time(std::chrono::time_point<std::chrono::steady_clock> chrono_clock_start,
                     std::chrono::time_point<std::chrono::steady_clock> chrono_clock_end){
    
    auto chromo_clock_diff = chrono_clock_end - chrono_clock_start;
    int secs;
    int hours;
    int mins;
    cout << "\nchrono::steady_clock runtime: \n";
    hours = std::chrono::duration_cast<std::chrono::hours>(chromo_clock_diff).count();
    cout << hours << "\tHours\n";
    mins  = std::chrono::duration_cast<std::chrono::minutes>(chromo_clock_diff).count();
    cout << mins - hours * 60 << "\tMinutes\n";
    secs  = std::chrono::duration_cast<std::chrono::seconds>(chromo_clock_diff).count();
    cout << secs - mins * 60 << "\tSeconds\n";
}

void print_system_time(std::chrono::time_point<std::chrono::system_clock> chrono_clock_start,
                       std::chrono::time_point<std::chrono::system_clock> chrono_clock_end){
    
    auto chromo_clock_diff = chrono_clock_end - chrono_clock_start;
    int secs;
    int hours;
    int mins;
    cout << "\nchrono::system_clock runtime: \n";
    hours = std::chrono::duration_cast<std::chrono::hours>(chromo_clock_diff).count();
    cout << hours << "\tHours\n";
    mins  = std::chrono::duration_cast<std::chrono::minutes>(chromo_clock_diff).count();
    cout << mins - hours * 60 << "\tMinutes\n";
    secs  = std::chrono::duration_cast<std::chrono::seconds>(chromo_clock_diff).count();
    cout << secs - mins * 60 << "\tSeconds\n";
}
