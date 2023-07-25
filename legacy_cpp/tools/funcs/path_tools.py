#!/usr/bin/env python3
# -*- coding: utf-8 -*

import os
import re
import glob

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
        if folder_name!="":
            if self.folder_name_in_path(folder_name):
                self.replace_folder_name(folder_name)
            else:
                self.path += folder_name + "/"


def get_file_paths(project_name):
    
    paths = glob.glob(f'{project_name}/**/', recursive=True)
    paths = [f for f in paths if '/buffs/' not in f]
    paths = [f for f in paths if 'time' in f]
    paths = [f+f.split('/')[-2] for f in paths]
    paths.sort()
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
