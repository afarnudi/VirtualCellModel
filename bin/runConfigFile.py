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
        self.bending = args.bending
        self.spring = args.spring
        self.platformID = args.platformID
        self.platformDeviceID = args.platformDeviceID
        self.rescale = 1 / args.rescale_elastic_properties

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

    def get_bending_coefficient(self):
        """retun bending coefficient"""
        return self.bending * self.rescale

    def get_spring_coefficient(self):
        """retun spring_coefficient"""
        return self.spring * self.rescale


def get_float_string(coef):
    """retun spring_coefficient"""
    coef = "{:.3f}".format(coef)
    while coef.endswith("0"):
        coef = coef[:-1]
    if coef.endswith("."):
        coef = coef[:-1]
    return coef


def configfile_gen_mem_ulms(seed, inputs, config_file_name):
    """modify excisting configfile"""
    print(seed)
    hostname = platform.node()
    edited_config_file_name = (
        config_file_name[:-4]
        + f"_{hostname}_{inputs.platformID}_{inputs.platformDeviceID}.txt"
    )

    temp = get_float_string(inputs.get_temperature())
    res = get_float_string(inputs.get_mesh_resolution())
    spring = get_float_string(inputs.get_spring_coefficient())
    bend = get_float_string(inputs.get_bending_coefficient())
    scale = get_float_string(inputs.get_scale())
    with open(config_file_name) as f:
        with open(edited_config_file_name, "w") as f1:
            for line in f:
                words = line.split(" ")
                if words[0] == "ProjectName":
                    line = line[:-1] + f"/{seed}\n"
                    line = re.sub(r"/\d+k", f"/{temp}k", line)
                    line = re.sub(r"/mesh_\d+", f"/mesh_{res}", line)
                    line = re.sub(r"/spring_\d+(\.\d*)?", f"/spring_{spring}", line)
                    line = re.sub(r"/bend_\d+(\.\d*)?", f"/bend_{bend}", line)
                    line = re.sub(r"/r_\d+", f"/r_{scale}", line)
                if words[0] == "MeshFile":
                    line = "MeshFile Mesh/icospheres/" + inputs.get_mesh_path()
                    line = line[:-5] + str(seed) + ".ply" + "\n"
                if words[0] == "Scale":
                    line = "Scale " + scale + "\n"
                if words[0] == "TemperatureinKelvin":
                    line = "TemperatureinKelvin " + temp + "\n"
                if words[0] == "SpringCoeff":
                    line = "SpringCoeff {:.4}\n".format(spring)
                if words[0] == "BendingCoeff":
                    line = "BendingCoeff {}\n".format(bend)
                f1.write(line)
    return edited_config_file_name


def get_arguments():
    """retunr arguments obtained with the argparse lib"""
    parser = argparse.ArgumentParser(
        description="Use VCM to run Memebrane samples with optional parameters."
    )
    parser.add_argument(
        "-p",
        "--platformID",
        help="Machine ID assigned to different platforms (OpenCL, CPU, etc)",
        default=2,
        type=int,
        metavar="\b",
    )
    parser.add_argument(
        "-c",
        "--config_file",
        help="config file name.",
        required=True,
        type=str,
        metavar="\b",
    )
    parser.add_argument(
        "--bending", help="Bending coefficient.", type=float, default=100, metavar="\b",
    )
    parser.add_argument(
        "--spring", help="Spring coefficient.", type=float, default=0.06, metavar="\b",
    )
    parser.add_argument(
        "-m",
        "--mesh_res",
        help="mesh resolution, M, with 10*M**2 + 2 nodes.",
        default=50,
        choices=[10, 30, 40, 50],
        type=int,
        metavar="\b",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        help="Simulation temperature",
        default=300,
        type=float,
        metavar="\b",
    )
    parser.add_argument(
        "-s",
        "--scale",
        help="scale of the membrane",
        default=2030,
        type=float,
        metavar="\b",
    )
    parser.add_argument(
        "-v",
        "--version",
        help="clang version of compiled VCM",
        default="clang++9",
        type=str,
        metavar="\b",
    )

    parser.add_argument(
        "-d",
        "--platformDeviceID",
        help="ID of device to use",
        default=0,
        type=int,
        metavar="\b",
    )
    parser.add_argument(
        "-r",
        "--rescale-elastic-properties",
        help="devide spring and bending rigidity by this number",
        default=1,
        type=float,
        metavar="\b",
    )
    parser.add_argument(
        "--seed_range",
        help="range of the seed. example: --seed_range 2 6, seeds 2,3,4, and 5 will be included.",
        required=True,
        type=int,
        nargs=2,
        metavar="\b",
    )
    return parser.parse_args()


def main():

    args = get_arguments()
    inputs = Inputs(args)

    config_file_name = args.config_file
    for seed in range(args.seed_range[0], args.seed_range[1]):
        edited_config_file_name = configfile_gen_mem_ulms(
            seed, inputs, config_file_name
        )
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
