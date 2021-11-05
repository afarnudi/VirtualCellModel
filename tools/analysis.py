#!/usr/bin/env python3
# -*- coding: utf-8 -*

import os
import argparse

import funcs.path_tools as pt



def main():
    parser = argparse.ArgumentParser(
        "Membrane analysis",
        description="Setup VCM.analysis to analyse membrane trajectories.",
    )
    parser.add_argument(
        "-t",
        "--tmux_max_session",
        help="Set max tmux session to use",
        type=int,
        default=4,
    )
    parser.add_argument(
        "-p",
        "--project_name",
        help="Set project directory.",
        type=str,
        default="MemTest",
    )
    parser.add_argument(
        "--clean", help="Find and delete all tmuxInProgress files", action="store_true"
    )
    args = parser.parse_args()

    project_name = os.getcwd() + "/"+ "Results/" + args.project_name
    tmux_session_count = 0
    tmux_max_sessions = args.tmux_max_session

    file_paths = pt.get_file_paths(project_name)
    os.system("tmux kill-server")
    for file_path in file_paths:
        ulm_file = file_path + "_mem0_Dulmts.txt"
        if not os.path.exists(ulm_file):
            tmux_session_in_progress = file_path + "_tmuxInProgress"
            if not args.clean:
                if not os.path.exists(tmux_session_in_progress):
                    if tmux_session_count < tmux_max_sessions:

                        seed = pt.get_seed_from_path(file_path)
                        mesh_path = pt.get_mesh_path(file_path, seed)
                        analysis_command = f"./VCM_AnalysisClang9 --analysis 3 --lmax 70 --framelimits 3,0 --meshfile {mesh_path} --ext _Dulmts.txt --filepath {file_path}.xyz"
                        try:
                            os.system(f"tmux new -s {tmux_session_count} -d")
                        except:
                            print(f"tmux session {tmux_session_count} already exists")
                        os.system(
                            f'tmux send -t {tmux_session_count} "{analysis_command}" C-m'
                        )
                        os.system(f"touch {tmux_session_in_progress}")
                        tmux_session_count += 1
            else:
                if os.path.exists(tmux_session_in_progress):
                    os.system(f"rm {tmux_session_in_progress}")

    print(f"{tmux_session_count} (out of {tmux_max_sessions}) sessions created")


if __name__ == "__main__":
    main()
