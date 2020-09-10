//
//  main.cpp
//  NodesOnSphere
//
//  Created by Ali Farnudi on 09/09/2020.
//  Copyright © 2020 Ali Farnudi. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <algorithm>

#include "Openmm.h"
#include "Openmm_structs.h"

using namespace std;

void generate_random_coordinates(int         N,
                                 MyAtomInfo *atoms);
void myWritePDBFrame(int frameNum,
                     double timeInPs,
                     const MyAtomInfo atoms[],
                     string traj_name);
void myGetOpenMMState(MyOpenMMData* omm,
                      double& timeInPs,
                      MyAtomInfo atoms[]);
MyOpenMMData* init_openmm(MyAtomInfo *atoms,
                          double stepSizeInFs);
void export_xyz(MyAtomInfo* atoms,
                string traj_name);
int main(int argc, const char * argv[]) {
    
    int N = 100, EndOfList=-1;
    string pdbname = "NodesOnSphere.pdb";
    
    
    if (argc>1) {
        N = stoi(argv[1]);
    }
    //extra nodes at the origin to keep everything constrained;
    N+=2;
    MyAtomInfo* atoms     = new MyAtomInfo[N+1];
    atoms[N].type         =EndOfList;
    generate_random_coordinates(N, atoms);
    
    
    float progressp=0;
    double Step_Size_In_Fs =1;
    int savetime = 0;
    double Report_Interval_In_Fs= 10000;
    double Simulation_Time_In_Ps =1000;
    int NumSilentSteps = (int)(Report_Interval_In_Fs / Step_Size_In_Fs + 0.5);
    try {
        MyOpenMMData* omm = new MyOpenMMData();
        omm = init_openmm(atoms, Step_Size_In_Fs);
        string platformName = omm->context->getPlatform().getName();
        
        cout<<"REMARK  Using OpenMM platform ";
        cout<<platformName.c_str()<<endl;
        myWritePDBFrame(0, 0, atoms, pdbname);
        
        for (int frame=1; ; ++frame) {
            
            double time;
            
            myGetOpenMMState(omm, time, atoms);
            
            if (int(time*1000/Step_Size_In_Fs) > savetime  ) {
                myWritePDBFrame(frame, time, atoms, pdbname);
                savetime += NumSilentSteps;
            }
            if (time >= Simulation_Time_In_Ps){
                break;
            }
            if ( omm->integrator != NULL ) {
                
                omm->integrator->step(NumSilentSteps);
                //        omm->context->computeVirtualSites();
            } else {
                
                omm->Lintegrator->step(NumSilentSteps);
                //        omm->context->computeVirtualSites();
            }
            
            if (100*time/Simulation_Time_In_Ps>progressp){
                printf("[ %2.1f ] time: %4.1f Ps [out of %4.1f Ps]    \r",100*time/Simulation_Time_In_Ps, time, Simulation_Time_In_Ps);
                cout<< std::flush;
                progressp =  int(1000*time/Simulation_Time_In_Ps)/10. + 0.1;
                
            }
            
        }
        
        cout<<"[ 100% ]\t time: "<<Simulation_Time_In_Ps<<"Ps\n";
            
        
        
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 0;
    }
    
    export_xyz(atoms,pdbname);
    
    return 0;
}

void export_xyz(MyAtomInfo* atoms,
                string traj_name){
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"w");
    int index=0;
    for (int n=1; atoms[n].type != EndOfList; ++n){
        
        fprintf(pFile,"%8.3f %8.3f %8.3f\n",
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2]);
    }
    fclose (pFile);
}

using OpenMM::Vec3;
MyOpenMMData* init_openmm(MyAtomInfo *atoms,
                          double stepSizeInFs){
    
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
    MyOpenMMData*       omm = new MyOpenMMData();
    OpenMM::System&     system = *(omm->system = new OpenMM::System());
    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
    string potential = "0.001*(sigma/r)^6";
    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
    ExcludedVolumes[0]->addGlobalParameter("sigma",   2*atoms[0].radius );
    ExcludedVolumes[0]->setCutoffDistance( 1.5 );
    system.addForce(ExcludedVolumes[0]);
    
    int EndOfList=-1;
    
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    
    for (int n=3; atoms[n].type != EndOfList; ++n) {
        HarmonicBond->addBond(0, n, 1, 20000);
    }
    system.addForce(HarmonicBond);
    
    
    std::vector<Vec3> initialPosInNm;
    std::vector<Vec3> initialVelInNmperPs;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        //        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atoms[n].mass);
        
        const Vec3 posInNm(atoms[n].initPosInNm[0],
                           atoms[n].initPosInNm[1],
                           atoms[n].initPosInNm[2]);
        const Vec3 velocityInNmperPs(atoms[n].velocityInNmperPs[0],
                                     atoms[n].velocityInNmperPs[1],
                                     atoms[n].velocityInNmperPs[2]);
        
        
        initialPosInNm.push_back(posInNm);
        initialVelInNmperPs.push_back(velocityInNmperPs);
        
        //add particles to the excluded volume force. The number of particles should be equal to the number particles in the system. The exluded interaction lists should be defined afterwards.
        ExcludedVolumes[0]->addParticle();
    }
    omm->EV = ExcludedVolumes;
    omm->harmonic = HarmonicBond;
    
    
    
    //cout<<"platform default directory path = "<<OpenMM::Platform::getDefaultPluginsDirectory()<<endl;
    //Listing the names of all available platforms.
    cout<<"\nOpenMM available platforms:\n"<<"Index Name \t  Speed (Estimated)\n";
    for (int i = 0; i < OpenMM::Platform::getNumPlatforms(); i++) {
        OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
        cout<<" ("<<std::to_string(i)<<")  "<<platform.getName().c_str()<<
        "\t   "<<platform.getSpeed()<<endl;
    }
    int platform_id=0;
    cout<<"Please choose a pltform (index): \n";
    std::cin>>platform_id;
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platform_id);
    
    
    std::vector<std::map<std::string, std::string> > device_properties;
    if (platform.getName() == "OpenCL") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["OpenCLPlatformIndex"]=std::to_string(i);
                    temp_device_properties["OpenCLDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<counter<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "OpenCLPlatformIndex") {
                            continue;
                        } else {
                            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                        }
                    }
                    cout<<"------------------------"<<endl;
                    counter++;
                    device_properties.push_back(temp_device_properties);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CUDA") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["CudaPlatformIndex"]=std::to_string(i);
                    temp_device_properties["CudaDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<counter<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "CUDAPlatformIndex") {
                            continue;
                        } else {
                            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                        }
                    }
                    cout<<"------------------------"<<endl;
                    counter++;
                    device_properties.push_back(temp_device_properties);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CPU") {
        OpenMM::System temp_system;
        temp_system.addParticle(1.0);
        OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
        OpenMM::Context temp_context(temp_system, temp_inegrator, platform);
        std::vector<std::string> platform_devices = platform.getPropertyNames();
        cout<<"CPU"<<" properties:\n";
        for (auto & name : platform_devices){
            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
        }
        cout<<"------------------------"<<endl;
    }
    
    
    int device_id=0;
    if (device_properties.size()>1) {
        cout<<"Please choose a device (index): \n";
        std::cin>>device_id;
    }
//    omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    omm->Lintegrator = new OpenMM::LangevinIntegrator(0,
                                                      0.01,
                                                      stepSizeInFs * OpenMM::PsPerFs);
    
    if (platform.getName() == "CPU" || platform_id==0) {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform);
        }
        
        
        
    } else {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, device_properties[device_id]);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform, device_properties[device_id]);
        }
    }
    
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    
    return omm;
}


void generate_random_coordinates(int N, MyAtomInfo* atoms){
    
    for (int i=0; i<3; i++) {
        atoms[i].type=0;
        atoms[i].pdb = "memb";
        atoms[i].energyInKJ = 0;
        atoms[i].symbol = 'M';
        atoms[i].velocityInNmperPs[0]=0;
        atoms[i].velocityInNmperPs[1]=0;
        atoms[i].velocityInNmperPs[2]=0;
        atoms[i].radius=sqrt(4./N);
    }
    atoms[1].initPosInNm[0]=0;
    atoms[1].initPosInNm[1]=0;
    atoms[1].initPosInNm[2]=0;
    atoms[1].posInNm[0]=0;
    atoms[1].posInNm[1]=0;
    atoms[1].posInNm[2]=0;
    atoms[1].mass=0;
    
    atoms[1].initPosInNm[0]=0;
    atoms[1].initPosInNm[1]=0;
    atoms[1].initPosInNm[2]=-1;
    atoms[1].posInNm[0]=0;
    atoms[1].posInNm[1]=0;
    atoms[1].posInNm[2]=-1;
    atoms[1].mass=0;
    
    atoms[2].initPosInNm[0]=0;
    atoms[2].initPosInNm[1]=0;
    atoms[2].initPosInNm[2]=1;
    atoms[2].posInNm[0]=0;
    atoms[2].posInNm[1]=0;
    atoms[2].posInNm[2]=1;
    
    for (int i=3; i<N; i++) {
        vector<double> xyz(3,0);
        double theta=((double)rand()/(double)RAND_MAX)*M_PI;
        double phi=((double)rand()/(double)RAND_MAX)*2*M_PI;
        
        double x = sin(theta)*cos(phi);
        double y = sin(theta)*sin(phi);
        double z = cos(theta);
        double dist=0;
        bool good = true;
//        cout<<atoms[0].radius<<endl;
        for (int j=1; j<i; j++) {
            dist= sqrt( (x-atoms[j].posInNm[0])*(x-atoms[j].posInNm[0]) + (y-atoms[j].posInNm[1])*(y-atoms[j].posInNm[1]) + (z-atoms[j].posInNm[2])*(z-atoms[j].posInNm[2]) );
//            cout<<dist<<endl;
            if (dist<atoms[0].radius) {
                good=false;
            }
        }
        if (!good) {
            i--;
            continue;
        }
        atoms[i].type=0;
        atoms[i].pdb = "memb";
        atoms[i].energyInKJ = 0;
        atoms[i].symbol = 'M';
        atoms[i].initPosInNm[0]=sin(theta)*cos(phi);
        atoms[i].initPosInNm[1]=sin(theta)*sin(phi);
        atoms[i].initPosInNm[2]=cos(theta);
        atoms[i].posInNm[0]=sin(theta)*cos(phi);
        atoms[i].posInNm[1]=sin(theta)*sin(phi);
        atoms[i].posInNm[2]=cos(theta);
        atoms[i].velocityInNmperPs[0]=0;
        atoms[i].velocityInNmperPs[1]=0;
        atoms[i].velocityInNmperPs[2]=0;
        atoms[i].mass=1;
        atoms[i].radius=sqrt(4./N);
    }
    
    
    
}

void myWritePDBFrame(int frameNum,
                     double timeInPs,
                     const MyAtomInfo atoms[],
                     std::string traj_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps\n",
            timeInPs);
    int index=0;
    
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    //    double occ=1;
    for (int n=1; atoms[n].type != EndOfList; ++n){
        
        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                atoms[n].pdb,
                chain[index],
                double(index),
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2],
                atoms[n].stretching_energy,
                atoms[n].energyInKJ);
    }
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}

void myGetOpenMMState(MyOpenMMData* omm,
                      double& timeInPs,
                      MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    
    const OpenMM::State state = omm->context->getState(infoMask,false);
    
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInNm[j] = positionsInNm[i][j];
        }
    }
}

