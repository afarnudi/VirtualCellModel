#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <algorithm>
#include <chrono>
#include <fstream>

#include "OpenMM.h"
#include "OpenMM_structs.h"

using namespace std;

void generate_random_coordinates(int N, MyAtomInfo *atoms, int seed);
void myWritePDBFrame(int frameNum, double timeInPs, const MyAtomInfo atoms[], string traj_name);
void writeXYZFrame  (int atom_count,
                     const MyAtomInfo atoms[],
                     double             time,
                     double             energyInKJ,
                     double             potential_energyInKJ,
                     std::string        traj_name
                     );
void myGetOpenMMState(MyOpenMMData* omm, double& timeInPs, MyAtomInfo atoms[]);
MyOpenMMData* init_openmm(MyAtomInfo *atoms, double stepSizeInFs, int seed, int &platformID, int &platformDeviceID);
void export_xyz(MyAtomInfo* atoms, string traj_name);
void print_wall_clock_time(double sim_duration_per_sec);
void print_real_time(chrono::time_point<chrono::steady_clock> chrono_clock_start,
                     chrono::time_point<chrono::steady_clock> chrono_clock_end);
void print_system_time(chrono::time_point<chrono::system_clock> chrono_clock_start,
                       chrono::time_point<chrono::system_clock> chrono_clock_end);

PlatformInfo get_platform_info(int &platformID, int &platformDeviceID);




int main(int argc, const char * argv[]) {
    
    auto chrono_clock_start = chrono::steady_clock::now();
    auto chrono_sys_clock_start = chrono::system_clock::now();
    clock_t tStart = clock();//Time the programme
    
    int N = 1002, EndOfList=-1;
    int seed=0;
    int platformID = -1;
    int platformDeviceID =-1;


    if (argc>1) {
        N = stoi(argv[1]);
        seed = stoi(argv[2]);
    }
    if (argc>3){
    	platformID = stoi(argv[3]);
	platformDeviceID = stoi(argv[4]); 
    }
    //cout<<N<<" "<<seed<<" "<<platformID<<" "<<platformDeviceID<<endl;
    string trajname = "USphere_"+to_string(N)+"d_s"+to_string(seed);
    string pdbname = trajname+".pdb";
    string xyzname = trajname+".xyz";
    //extra node at the origin to keep everything constrained;
    N++;
    //extra -1 at the end of list
    MyAtomInfo* atoms     = new MyAtomInfo[N+1];
    atoms[N].type         =EndOfList;
    generate_random_coordinates(N, atoms, seed);
    
    
    
    float progressp=0;
    double Step_Size_In_Fs =1;
    int savetime = 0;
    double Simulation_Time_In_Ps =200;
    double Report_Interval_In_Fs= 10000;
    int NumSilentSteps = (int)(Report_Interval_In_Fs / Step_Size_In_Fs + 0.5);
    try {
        MyOpenMMData* omm = new MyOpenMMData();
        omm = init_openmm(atoms, Step_Size_In_Fs, seed, platformID, platformDeviceID);
        string platformName = omm->context->getPlatform().getName();
        
        
        
        tStart = clock();
        chrono_clock_start = chrono::steady_clock::now();
        chrono_sys_clock_start = chrono::system_clock::now();
        
        cout<<"REMARK  Using OpenMM platform ";
        cout<<platformName.c_str()<<endl;
//        myWritePDBFrame(0, 0, atoms, pdbname);
        writeXYZFrame  (N, atoms, 0, 0, 0, xyzname);
        
        double minimiseTolerance =1;
        int minimiseMaxIterations=0;
        
        cout<<"Minimising the coordinates before running the simulation."<<endl;
        cout<<"Minimisation parameters:"<<endl;
        cout<<"\tMinimisation Tolerance: "<<minimiseTolerance<<endl;
        cout<<"\tMinimisation Max Iterations: "<<minimiseMaxIterations<<endl<<endl;
        OpenMM::LocalEnergyMinimizer::minimize(*(omm->context), minimiseTolerance, minimiseMaxIterations);
        double fake_time, fake_energyInKJ, fake_potential_energyInKJ;
        myGetOpenMMState(omm, fake_time, atoms);
        writeXYZFrame  (N, atoms, 0, 0, 0, xyzname);
        cout<<"\tMinimisation finished."<<endl<<endl;


        for (int frame=2; ; ++frame) {

            double time;

            myGetOpenMMState(omm, time, atoms);

            if (int(time*1000/Step_Size_In_Fs) > savetime  ) {
//                myWritePDBFrame(frame, time, atoms, pdbname);
                writeXYZFrame  (N, atoms, time, 0, 0, xyzname);
                savetime += NumSilentSteps;
            }
            if (time >= Simulation_Time_In_Ps){
                break;
            }
            if ( omm->integrator != NULL ) {

                omm->integrator->step(NumSilentSteps);
                //        omm->context->computeVirtualSites();
            } else {

                omm->LangevinIntegrator->step(NumSilentSteps);
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

        print_wall_clock_time((double)((clock() - tStart)/CLOCKS_PER_SEC));
        print_real_time(chrono_clock_start, chrono::steady_clock::now());
        print_system_time(chrono_sys_clock_start, chrono::system_clock::now());
        
//    export_xyz(atoms,pdbname);
    
    return 0;
}


void export_xyz(MyAtomInfo* atoms,
                string traj_name){
    traj_name.erase(traj_name.end()-3,traj_name.end());
    traj_name+="xyz";
    //traj_name = "NodeGenMems/" + traj_name;
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"w");
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
                          double stepSizeInFs,
                          int seed,
			  int &platformID,
			  int &platformDeviceID){
    const string cbp_plugin_location="/scratch/alifarnudi/local/openmm/lib/plugins";
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
//    OpenMM::Platform::loadPluginsFromDirectory(cbp_plugin_location);
    MyOpenMMData*       omm = new MyOpenMMData();
    OpenMM::System&     system = *(omm->system = new OpenMM::System());
    
    
    int EndOfList=-1;
    int num_of_atoms=0;
    double k_spring=200;
    
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    
    for (int n=1; atoms[n].type != EndOfList; ++n) {
        num_of_atoms++;
        HarmonicBond->addBond(0, n, 1, k_spring);
    }
    system.addForce(HarmonicBond);
    
    
//    cout<<num_of_atoms<<endl;
    
    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
//    string potential = to_string(k_spring*4*M_PI/(100*num_of_atoms))+"*(1/r)^2";
    double sigma = 2*sqrt((4.*M_PI)/num_of_atoms);
//    string sigma =to_string(sigma);
    string potential = to_string(4*k_spring/10000000.)+"*(("+to_string(sigma)+"/r)^8-("+to_string(sigma)+"/r)^4+0.25)";
    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
//
//    string potential = "0.001*(sigma/r)^6";
//    ExcludedVolumes.push_back(new OpenMM::CustomNonbondedForce(potential));
//    ExcludedVolumes[0]->addGlobalParameter("sigma",   2*atoms[0].radius );
    
    ExcludedVolumes[0]->setCutoffDistance( sigma*pow(2,1./4.) );
    system.addForce(ExcludedVolumes[0]);
    
    
    
    
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
    
    PlatformInfo platforminfo = get_platform_info(platformID, platformDeviceID);
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
        
    //    omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    omm->LangevinIntegrator = new OpenMM::LangevinIntegrator(0,
                                                             0.1,
                                                             stepSizeInFs * OpenMM::PsPerFs);
    omm->LangevinIntegrator->setRandomNumberSeed(seed);
    
    if (platform.getName() == "CPU" || platformID==0) {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinIntegrator, platform);
        }
        
        
        
    } else {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinIntegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        }
    }
//    exit(0);
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    
    return omm;
}


void generate_random_coordinates(int N, MyAtomInfo* atoms, int seed){
    srand(seed);
    char pdblabel[]={'m','e','m','b'};
    
    atoms[0].type=0;
    atoms[0].pdb = pdblabel;
    atoms[0].energyInKJ = 0;
    atoms[0].symbol = 'M';
    atoms[0].velocityInNmperPs[0]=0;
    atoms[0].velocityInNmperPs[1]=0;
    atoms[0].velocityInNmperPs[2]=0;
    atoms[0].radius=sqrt(4./N);
    atoms[0].initPosInNm[0]=0;
    atoms[0].initPosInNm[1]=0;
    atoms[0].initPosInNm[2]=0;
    atoms[0].posInNm[0]=0;
    atoms[0].posInNm[1]=0;
    atoms[0].posInNm[2]=0;
    atoms[0].mass=0;
    
    for (int i=1; i<N; i++) {
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
        atoms[i].pdb = pdblabel;
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
    atoms[N].type=-1;
    
    
    
}

void writeXYZFrame  (int atom_count,
                     const MyAtomInfo atoms[],
                     double             time,
                     double             energyInKJ,
                     double             potential_energyInKJ,
                     std::string        traj_name
                     )
{
//    traj_name= generalParameters.buffer_file_name+".xyz_buff";
    ofstream writexyz(traj_name.c_str(), ios_base::app);
    writexyz<<atom_count<<endl;
    writexyz<<"timePs ";
    writexyz<<time<<setprecision(16);
    writexyz<<" potential_energy_inKJpermol ";
    writexyz<<potential_energyInKJ<<setprecision(16);
    writexyz<<" energy_inKJpermol ";
    writexyz<<energyInKJ<<setprecision(16)<<endl;
    for (int n=0; atoms[n].type != -1; n++) {
        writexyz<<"Mem0\t"<<atoms[n].posInNm[0]<<"\t"<<atoms[n].posInNm[1]<<"\t"<<atoms[n].posInNm[2]<<setprecision(16)<<"\n";
    }
    writexyz.close();
    
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

PlatformInfo get_platform_info(int &platformID, int &platformDeviceID)
{
    int stepSizeInFs =1;
    int list_depth=20;
    PlatformInfo platforminfo;
    
    
    //cout<<"platform default directory path = "<<OpenMM::Platform::getDefaultPluginsDirectory()<<endl;
    //Listing the names of all available platforms.
//    cout<<"Here";
    cout<<"\nOpenMM available platforms:\n"<<"Index Name \t  Speed (Estimated)\n";
    for (int i = 0; i < OpenMM::Platform::getNumPlatforms(); i++) {
        OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
        cout<<" ("<<std::to_string(i)<<")  "<<platform.getName().c_str()<<
        "\t   "<<platform.getSpeed()<<endl;
    }
    platforminfo.platform_id=0;
    if (platformID!=-1) {
        platforminfo.platform_id = platformID;
    } else {
        cout<<"Please choose a pltform (index): \n";
        std::cin>>platforminfo.platform_id;
    }
    
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
    
    //    std::vector<std::string> device_properties_report;
    //    std::vector<std::map<std::string, std::string> > device_properties;
    if (platform.getName() == "OpenCL") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int i=0; i<list_depth; i++) {
            for (int j=0; j<list_depth; j++) {
                vector<string> data={"single","double", "mixed"};
                for (vector<string>::iterator precision=data.begin(); precision!=data.end(); ++precision) {
                    try {
                        std::string report;
                        std::map<std::string, std::string> temp_device_properties;
                        temp_device_properties["OpenCLPlatformIndex"]=std::to_string(i);
                        temp_device_properties["OpenCLDeviceIndex"]=std::to_string(j);
                        temp_device_properties["Precision"]=*precision;
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
                                report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                                cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                            }
                        }
                        cout<<"------------------------"<<endl;
                        counter++;
                        platforminfo.device_properties.push_back(temp_device_properties);
                        platforminfo.device_properties_report.push_back(report);
                    } catch (const std::exception& e) {
                        
                    }
                }
            }
        }
    } else if (platform.getName() == "CUDA") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int j=0; j<2*list_depth; j++) {
            vector<string> data={"single","double", "mixed"};
            for (vector<string>::iterator precision=data.begin(); precision!=data.end(); ++precision) {
                try {
                    std::string report;
                    std::map<std::string, std::string> temp_device_properties;
                    //                    temp_device_properties["CudaPlatformIndex"]=std::to_string(i);
                    temp_device_properties["DeviceIndex"]=std::to_string(j);
                    temp_device_properties["CudaPrecision"]=*precision;
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<counter<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex") {
                            continue;
                        } else {
                            report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                        }
                    }
                    cout<<"------------------------"<<endl;
                    counter++;
                    platforminfo.device_properties.push_back(temp_device_properties);
                    platforminfo.device_properties_report.push_back(report);
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
    
    platforminfo.platform_device_id=0;
    if (platforminfo.device_properties.size()>1) {
        if (platformDeviceID!=-1) {
            platforminfo.platform_device_id = platformDeviceID;
        } else {
            cout<<"Please choose a device (index): \n";
            std::cin>>platforminfo.platform_device_id;
        }
        
    }
    
    return platforminfo;
}
