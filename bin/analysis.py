#!/usr/bin/env python3
# -*- coding: utf-8 -*

import os
import re
import argparse

def get_folder_names_and_values(path, prefix, suffix, data_type):
    folders = os.listdir(path)
    patern = prefix + r"[-+]?\d*\.\d+" + suffix + "|" + prefix + r"\d+" + suffix
    names = [f for f in folders if re.search(patern, f)]
    values = [data_type(re.findall(r"[-+]?\d*\.\d+|\d+", f)[0]) for f in names]
    return [(a, b) for a, b in zip(names, values)]


class Path:
    def __init__(self, directory):
        if directory[-1] != "/":
            directory += "/"
        self.path = directory

    def get_dir(self,):
        return self.path

    def filter_empty_strings(list_of_strings):
        return [f for f in list_of_strings if f != ""]

    def get_value_and_index(list_of_strings):
        value_index = None
        for part in list_of_strings:
            try:
                value = float(part)
                try:
                    value_index = list_of_strings.index(str(value))
                except ValueError:
                    value = int(part)
                    value_index = list_of_strings.index(str(value))
            except:
                # I am looking for the number part
                pass
        return value, value_index

    def get_keywords(list_of_strings):
        parts = re.split(r"(\d+(?:\.\d+)?)", list_of_strings)
        result = []
        for string in parts:
            split = string.split("/")
            split = Path.filter_empty_strings(split)
            for key in split:
                result.append(key)
        return result

    def folder_name_in_path(self, folder_name):
        parts = Path.get_keywords(folder_name)
        # print(parts)
        vlue, value_index = Path.get_value_and_index(parts)
        dir_parts = Path.get_keywords(self.path)
        # print(dir_parts)
        num_of_keys_in_path = 0
        for i in range(len(parts)):
            if i != value_index:
                num_of_keys_in_path += parts[i] in dir_parts
        return num_of_keys_in_path == len(parts) - 1

    def replace_folder_name(self, folder_name):
        parts = Path.get_keywords(folder_name)
        vlue, value_index = Path.get_value_and_index(parts)
        try:
            prefix = parts[:value_index][0]
        except IndexError:
            prefix = ""

        try:
            suffix = parts[value_index + 1 :][0]
        except IndexError:
            suffix = ""
        self.path = re.sub(
            r"{}[-+]?\d*\.\d+{}/|/{}\d+{}".format(prefix, suffix, prefix, suffix),
            "/" + folder_name + "/",
            self.path,
        )
        self.path = self.path.split("/" + folder_name + "/")[0]
        if self.path[-1] != "/":
            self.path += "/"
        self.path += folder_name + "/"

    def add_folder(self, folder_name):

        if self.folder_name_in_path(folder_name):
            self.replace_folder_name(folder_name)
        else:
            self.path += folder_name + "/"


def get_file_paths(project_name, results_path):
    paths = []
    project_path = Path(project_name)
    steps = get_folder_names_and_values(
        results_path + project_path.get_dir(), 
        "step_", 
        "", 
        int
    )
    for step in steps:
        step_name, step_value = step
        project_path.add_folder(step_name)
        frictions = get_folder_names_and_values(
            results_path + project_path.get_dir(), 
            "friction_", 
            "", 
            float
        )
        for friction in frictions:
            friction_name, friction_value = friction
            project_path.add_folder(friction_name)
            temperatures = get_folder_names_and_values(
                results_path + project_path.get_dir(), 
                "", 
                "k", 
                int
            )
            for temperature in temperatures:
                temperature_name, temperature_value = temperature
                project_path.add_folder(temperature_name)
                meshes = get_folder_names_and_values(
                    results_path + project_path.get_dir(), 
                    "mesh_", 
                    "", 
                    int
                )
                for mesh in meshes:
                    mesh_name, mesh_value = mesh
                    project_path.add_folder(mesh_name)
                    sigmas = get_folder_names_and_values(
                        results_path + project_path.get_dir(), 
                        "spring_", 
                        "", 
                        float
                    )
                    for sigma in sigmas:
                        sigma_name, sigma_value = sigma
                        project_path.add_folder(sigma_name)
                        kappas = get_folder_names_and_values(
                            results_path + project_path.get_dir(), 
                            "bend_", 
                            "", 
                            float
                        )
                        for kappa in kappas:
                            kappa_name, kappa_value = kappa
                            project_path.add_folder(kappa_name)
                            radii = get_folder_names_and_values(
                                results_path + project_path.get_dir(), 
                                "r_", 
                                "", 
                                float
                            )
                            for radius in radii:
                                radius_name, radius_value = radius
                                project_path.add_folder(radius_name)
                                samples = get_folder_names_and_values(
                                    results_path + project_path.get_dir(), 
                                    "", 
                                    "", 
                                    int
                                )
                                smps = []
                                for smp in samples:
                                    try:
                                        int(smp[0])
                                        smps.append(smp)
                                    except:
                                        #I don't need samples that are not integers
                                        continue
                                samples=smps
                                samples.sort()
                                for seed in samples:
                                    seed_name, seed_value = seed
                                    file_path = (
                                        os.getcwd()
                                        + "/"
                                        + results_path
                                        + project_path.get_dir()
                                        + f"{seed_name}/"
                                    )
                                    folders = os.listdir(file_path)
                                    folders.sort()
                                    paths.append(
                                        file_path + folders[-1] + "/" + folders[-1]
                                    )
    return paths


def get_param_value(path, prefix, suffix, data_type):
    patern = prefix + r"[-+]?\d*\.\d+" + suffix + "|" + prefix + r"\d+" + suffix
    param_string = re.findall(patern, path)[0]
    return data_type(re.findall(r"[-+]?\d*\.\d+|\d+", param_string)[0])


def get_mesh_path(path, seed):
    value = get_param_value(path, "mesh_", "", int)
    mesh_nodes = 10 * value ** 2 + 2
    return f"Mesh/icospheres/USphere_{mesh_nodes}d/USphere_{mesh_nodes}d_s{seed}.ply"


def get_seed_from_path(path):
    return int(path.split("/")[-3])


def main():
    parser = argparse.ArgumentParser(
        "Membrane analysis", description="Setup VCM.analysis to analyse membrane trajectories."
    )
    parser.add_argument("-t",
                        "--tmux_max_session", 
                        help="Set max tmux session to use",
                        type=int,
                        default=4)
    parser.add_argument("-p",
                        "--project_name", 
                        help="Set project directory.",
                        type=str,
                        default="MemTest")
    parser.add_argument("--clean", 
                        help="Find and delete all tmuxInProgress files", 
                        action= "store_true")
    args = parser.parse_args()
    
    
    project_name = args.project_name
    results_path = "Results/"
    tmux_session_count = 0
    tmux_max_sessions = args.tmux_max_session

    file_paths = get_file_paths(project_name, results_path)
    for file_path in file_paths:
        ulm_file = file_path + "_mem0_Dulmts.txt"
        if not os.path.exists(ulm_file):
            tmux_session_in_progress = file_path + "_tmuxInProgress"
            if not args.clean:
                if not os.path.exists(tmux_session_in_progress):
                    if tmux_session_count < tmux_max_sessions:
                        os.system(f"touch {tmux_session_in_progress}")
                        seed = get_seed_from_path(file_path)
                        mesh_path = get_mesh_path(file_path, seed)
                        analysis_command = f"./VCM_AnalysisClang9 --analysis 3 --lmax 70 --framelimits 3,0 --meshfile {mesh_path} --ext _Dulmts.txt --filepath {file_path}.xyz"
                        os.system(f"tmux new -s {tmux_session_count} -d")
                        os.system(
                            f'tmux send -t {tmux_session_count} "{analysis_command}" C-m'
                        )
                        tmux_session_count += 1
            else:
                if os.path.exists(tmux_session_in_progress):
                    os.system(f"rm {tmux_session_in_progress}")

    print(f"{tmux_session_count} (out of {tmux_max_sessions}) sessions created")
if __name__ == "__main__":
    main()
