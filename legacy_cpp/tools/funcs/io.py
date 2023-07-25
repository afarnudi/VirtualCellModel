#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:19:43 2021

@author: ali
"""
import re
import os
import numpy as np
import pandas as pd
from .path_tools import get_folder_names_and_values, get_value_from_path
from .path_tools import get_file_paths
from .ulm_analysis import process_spherical_harmonics_amplitudes
from .ulm_analysis import equatorial_mode_calculator


def extract_paths(paths_with_samples):
    paths = []
    for path in paths_with_samples:
        path_no_sample = re.split(r"/\d+/", path)[0] + "/"
        if not path_no_sample in paths:
            paths.append(path_no_sample)
    return paths


def load_data(project_name):

    columns = [
        "gamma",
        "kappa",
        "Young2D",
        "R",
        "T",
        "ul",
        "ul_error",
        "um2",
        "um2_error",
        "mesh",
        "step_size",
        "friction",
    ]
    Gamma_column = []
    kappa_column = []
    Young2D_column = []
    R_column = []
    T_column = []
    mesh_column = []
    avg_column = []
    std_column = []
    um2_column = []
    um2_std_column = []
    step_column = []
    friction_column = []

    paths = get_file_paths(project_name)
    paths = extract_paths(paths)
    for path in paths:
        step_value        = get_value_from_path(path, "step_", "", int)
        friction_value    = get_value_from_path(path, "friction_", "", float)
        temperature_value = get_value_from_path(path, "", "k", int)
        mesh_value        = get_value_from_path(path, "mesh_", "", int)
        sigma_value       = get_value_from_path(path, "spring_", "", float)
        kappa_value       = get_value_from_path(path, "bend_", "", float)
        radius_value      = get_value_from_path(path, "r_", "", float)

        samples = get_folder_names_and_values(path, "", "", int)
        samples.sort()
        kappa_column.append(kappa_value * np.sqrt(3) / 2)
        Young2D_column.append(sigma_value * 2 / np.sqrt(3))
        R_column.append(radius_value)
        T_column.append(temperature_value)
        Gamma_column.append(sigma_value * radius_value ** 2 / kappa_value)
        mesh_column.append(mesh_value)
        step_column.append(step_value)
        friction_column.append(friction_value)

        ulmMatrixAvg = []
        ulmAvg = []
        um2Avg = []

        for sample in samples:
            sample_name, sample_value = sample

            folders = os.listdir(path + sample_name)
            fname = [f for f in folders if "2021" in f]
            fname.sort()

            (
                ulm_squared_avg,
                ulm_squared_std,
                ell_max,
                equitorial_amplitudes_avg,
                equitorial_amplitudes_std,
                ulm_squared_matrix_avg,
                ulm_squared_matrix_std,
            ) = process_spherical_harmonics_amplitudes(
                path + sample_name, file_names=fname[-1:], extension="_mem0_Dulmts.txt",
            )

            ulmMatrixAvg.append(ulm_squared_matrix_avg)
            ulmAvg.append(ulm_squared_avg)
            um2Avg.append(equitorial_amplitudes_avg)
        avg_column.append(np.mean(ulmAvg, axis=0))
        std_column.append(np.std(ulmAvg, axis=0) / np.sqrt(len(samples) - 1))

        AVGulmMatrixAvg = np.mean(ulmMatrixAvg, axis=0)
        STDulmMatrixAvg = np.std(ulmMatrixAvg, axis=0) / np.sqrt(len(samples) - 1)
        if mesh_value == 10:
            mMax = 14
            um2_column.append(equatorial_mode_calculator(AVGulmMatrixAvg, mMax))
            um2_std_column.append(equatorial_mode_calculator(STDulmMatrixAvg, mMax))
        else:
            um2_column.append(equatorial_mode_calculator(AVGulmMatrixAvg))
            um2_std_column.append(equatorial_mode_calculator(STDulmMatrixAvg))

    data = {
        columns[0]: Gamma_column,
        columns[1]: kappa_column,
        columns[2]: Young2D_column,
        columns[3]: R_column,
        columns[4]: T_column,
        columns[5]: avg_column,
        columns[6]: std_column,
        columns[7]: um2_column,
        columns[8]: um2_std_column,
        columns[9]: mesh_column,
        columns[10]: step_column,
        columns[11]: friction_column,
    }
    df = pd.DataFrame(data, columns=columns)
    return df
