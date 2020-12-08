#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;

void print_platform_info(void){
    const string cbp_plugin_location="/scratch/alifarnudi/local/openmm/lib/plugins";
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
    //OpenMM::Platform::loadPluginsFromDirectory(cbp_plugin_location);
    PlatformInfo platforminfo = get_platform_info();
    cout<<"device flags:\n--platformID "<<platforminfo.platform_id<<" --platformDeviceID "<<platforminfo.platform_device_id<<endl;
    
}

PlatformInfo get_platform_info(void)
{
    int stepSizeInFs =1;
    
    PlatformInfo platforminfo;
    
    
    //cout<<"platform default directory path = "<<OpenMM::Platform::getDefaultPluginsDirectory()<<endl;
    //Listing the names of all available platforms.
    cout<<TOMM<<"\nOpenMM available platforms:\n"<<TGRAY<<"Index Name \t  Speed (Estimated)\n"<<TRESET;
    for (int i = 0; i < OpenMM::Platform::getNumPlatforms(); i++) {
        OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
        cout<<" ("<<TBOLD<<std::to_string(i)<<TRESET<<")  "<<platform.getName().c_str()<<
        "\t   "<<platform.getSpeed()<<endl;
    }
    platforminfo.platform_id=0;
    cout<<"Please choose a pltform (index): \n"<<TFILE;
    std::cin>>platforminfo.platform_id;
    cout<<TRESET;
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
    
//    std::vector<std::string> device_properties_report;
//    std::vector<std::map<std::string, std::string> > device_properties;
    if (platform.getName() == "OpenCL") {
        cout<<"Available devices on the "<<TOCL<<platform.getName()<<TRESET<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::string report;
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["OpenCLPlatformIndex"]=std::to_string(i);
                    temp_device_properties["OpenCLDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<TBOLD<<TOCL<<counter<<TRESET<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "OpenCLPlatformIndex") {
                            continue;
                        } else {
                            report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                            cout<<"\t"<<name<<"\t"<<TOCL<<platform.getPropertyValue(temp_context, name)<<TRESET<<endl;
                        }
                    }
                    cout<<TRESET<<"------------------------"<<endl;
                    counter++;
                    platforminfo.device_properties.push_back(temp_device_properties);
                    platforminfo.device_properties_report.push_back(report);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CUDA") {
        cout<<"Available devices on the "<<TCUD<<platform.getName()<<TRESET<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::string report;
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["CudaPlatformIndex"]=std::to_string(i);
                    temp_device_properties["CudaDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<TBOLD<<TCUD<<counter<<TRESET<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "CUDAPlatformIndex") {
                            continue;
                        } else {
                            report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                            cout<<"\t"<<name<<"\t"<<TCUD<<platform.getPropertyValue(temp_context, name)<<TRESET<<endl;
                        }
                    }
                    cout<<TRESET<<"------------------------"<<endl;
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
        cout<<TCPU<<"CPU"<<TRESET<<" properties:\n";
        for (auto & name : platform_devices){
            cout<<"\t"<<name<<"\t"<<TCPU<<platform.getPropertyValue(temp_context, name)<<TRESET<<endl;
        }
        cout<<TRESET<<"------------------------"<<endl;
    }
    
    int device_id=0;
    if (platforminfo.device_properties.size()>1) {
        cout<<"Please choose a device (index): \n"<<TFILE;
        std::cin>>device_id;
        cout<<TRESET;
    }
    
    return platforminfo;
}

void get_platform_info(PlatformInfo &platforminfo)
{
    int stepSizeInFs =1;
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
    
    if (platform.getName() == "OpenCL") {
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::string report;
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["OpenCLPlatformIndex"]=std::to_string(i);
                    temp_device_properties["OpenCLDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "OpenCLPlatformIndex") {
                            continue;
                        } else {
                            report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                        }
                    }
                    counter++;
                    platforminfo.device_properties.push_back(temp_device_properties);
                    platforminfo.device_properties_report.push_back(report);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CUDA") {
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::string report;
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["CudaPlatformIndex"]=std::to_string(i);
                    temp_device_properties["CudaDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "CUDAPlatformIndex") {
                            continue;
                        } else {
                            report+="\t"+name+"\t"+platform.getPropertyValue(temp_context, name)+"\n";
                        }
                    }
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
    }
}

void generateHardwareReport (PlatformInfo platforminfo){
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
    int stepSizeInFs =1;
    generalParameters.hardwareReport ="Running on the "+platform.getName()+" platform:\n";
    if (platform.getName() != "CPU") {
        generalParameters.hardwareReport+=platforminfo.device_properties_report[platforminfo.platform_device_id]+"\n";
    } else {
        OpenMM::System tsystem;
        tsystem.addParticle(1.0);
        OpenMM::VerletIntegrator tinegrator(stepSizeInFs * OpenMM::PsPerFs);
        OpenMM::Context tcontext(tsystem, tinegrator, platform);
        std::vector<std::string> tplatform_devices = platform.getPropertyNames();
        for (auto & name : tplatform_devices){
            generalParameters.hardwareReport+="\t"+name+"\t"+platform.getPropertyValue(tcontext, name)+"\n";
        }
    }
    generalParameters.hardwareReport+="------------------------\n\n";
}
