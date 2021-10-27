import argparse
import sys
import time
import platform
import re
import os


class Inputs:
    """Parse argparse arguments"""

    def __init__(self, args):
        self.mesh_resolution = args.mesh_res
        self.mesh_path = Inputs.parse_mesh_resolution(args)
        self.temperature = Inputs.parse_temperature(args)
        self.scale = Inputs.parse_scale(args)

    def generate_mesh_path_from_resolution(x):
        return "USphere_{0}2d/USphere_{0}2d_s0.ply".format(x ** 2)

    def parse_mesh_resolution(args):
        """Parse resolution input to mesh path"""
        if args.mesh_res in [50, 10]:
            return Inputs.generate_mesh_path_from_resolution(args.mesh_res)
        else:
            raise Exception("only available options for mesh resolutions are 10 and 50")

    def parse_temperature(args):
        """assign temperature"""
        return float(args.temperature)

    def parse_scale(args):
        """assign scale"""
        return float(args.scale)

    def get_scale(self,):
        """return scale"""
        return float(self.scale)

    def get_temperature(self):
        """retun temperature in string format"""
        return int(self.temperature)

    def get_mesh_resolution(self):
        """retun mesh resolution (int)"""
        return self.mesh_resolution

    def get_mesh_path(self):
        """retun mesh path (str)"""
        return self.mesh_path


def configfile_gen_mem_ulms(index, inputs, config_file_name):
    """modify excisting configfile"""
    print(index)
    hostname = platform.node()
    edited_config_file_name = config_file_name[:-4] + f"_{hostname}.txt"
    with open(config_file_name) as f:
        with open(edited_config_file_name, "w") as f1:
            for line in f:
                words = line.split(" ")
                if words[0] == "ProjectName":
                    line = line[:-1] + f"/{index}\n"
                    line = re.sub(r"/\d+k", f"/{inputs.get_temperature()}k", line)
                    line = re.sub(
                        r"/mesh_\d+", f"/mesh_{inputs.get_mesh_resolution()}", line
                    )
                    line = re.sub(r"/r_\d+", f"/r_{inputs.get_scale()}", line)
                if words[0] == "MeshFile":
                    line = "MeshFile Mesh/icospheres/" + inputs.get_mesh_path()
                    line = line[:-5] + str(index) + ".ply" + "\n"
                if words[0] == "Scale":
                    line = "Scale " + str(inputs.get_scale()) + "\n"
                if words[0] == "TemperatureinKelvin":
                    line = "TemperatureinKelvin " + str(inputs.get_temperature()) + "\n"
                f1.write(line)
    return edited_config_file_name


def get_arguments():
    """retunr arguments obtained with the argparse lib"""
    parser = argparse.ArgumentParser(
        description="Use VCM to run Memebrane samples with optional parameters."
    )
    parser.add_argument(
        "-m", "--mesh-res",
        help="mesh resolution, 10 for 1002, and 50 for 25002 nodes.",
        default=50,
        choices=[10,50],
        type=int,
        metavar='',
    )
    parser.add_argument(
        "-t", "--temperature", 
        help="Simulation temperature", 
        default=300, 
        type=float,
        metavar='',
    )
    parser.add_argument(
        "-s", "--scale", 
        help="scale of the membrane", 
        default=100, 
        type=float,
        metavar='',
    )
    parser.add_argument(
        "-v",
        "--version",
        help="clang version of compiled VCM",
        default="clang++9",
        type=str,
        metavar='',
    )
    parser.add_argument(
        "-p",
        "--platformID",
        help="Machine ID assigned to different platforms, i.e. OpenCL, CPU, Cuda, ...",
        default=2,
        type=int,
        metavar='',
    )
    parser.add_argument(
        "-d", "--platformDeviceID", 
        help="ID of device to use", 
        default=0, 
        type=int,
        metavar='',
    )
    return parser.parse_args()


def main():

    args = get_arguments()
    inputs = Inputs(args)

    config_file_name = "membraneTest.txt"
    hostname = platform.node()
    for i in range(10):
        edited_config_file_name = configfile_gen_mem_ulms(i, inputs, config_file_name)
        if args.version == "clang++":
            os.system(
                f"./VCM -c {edited_config_file_name} --platformID {args.platformID} --platformDeviceID {args.platformDeviceID}"
            )
        elif args.version == "clang++9":
            os.system(
                f"./VCMclang9 -c {edited_config_file_name} --platformID {args.platformID} --platformDeviceID {args.platformDeviceID}"
            )

    print("Finished")


if __name__ == "__main__":
    main()
