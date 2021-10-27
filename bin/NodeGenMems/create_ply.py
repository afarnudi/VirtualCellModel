#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import argparse


def print_properties(tris, bonds, xyzs):
    print("number of nodes ", xyzs.shape[0])
    print("number of node pairs ", bonds.shape[0])
    print("number of triangles: ", tris.shape[0])
    degree = np.zeros(xyzs.shape[0], dtype=int)
    for i in bonds.reshape(bonds.shape[0] * bonds.shape[1]):
        degree[i - 1] += 1

    hist = {}
    degreenp = np.asarray(list(set(degree)))
    for i in degreenp:
        hist[i] = 0
    for i in degree:
        hist[i] += 1
    print("node degree:count  ", hist)
    print("avg degree:   ", np.mean(degreenp))


def get_bond_set(bonds):
    bond_set = {tuple(bonds[0])}
    for i in bonds:
        i.sort()
        bond_set.add(tuple(i))
    return bond_set


def get_triangle_list(bond_set):
    triangle_set = set()
    for i in bond_set:
        node1 = i[0]
        node2 = i[1]
        for j in bond_set:
            if node1 == j[0]:
                node3 = j[1]
                if node3 != node2:
                    lnode = [node2, node3]
                    lnode.sort()
                    if tuple(lnode) in bond_set:
                        ltri = [node1, node2, node3]
                        ltri.sort()
                        triangle_set.add(tuple(ltri))
            if node1 == j[1]:
                node3 = j[0]
                if node3 != node2:
                    lnode = [node2, node3]
                    lnode.sort()
                    if tuple(lnode) in bond_set:
                        ltri = [node1, node2, node3]
                        ltri.sort()
                        triangle_set.add(tuple(ltri))
    triangle_list = list(triangle_set)
    triangle_list.sort()
    return triangle_list


def get_bond_array(sphere_delaunay_bonds_output):
    bonds = np.loadtxt(sphere_delaunay_bonds_output, dtype=int)
    N = bonds.shape[0] // 3
    bonds = bonds.reshape((N, 3))
    bonds = bonds[:, :2]
    bonds[0].sort()
    return bonds


def build_triangle_from_bonds(sphere_delaunay_bonds_output):
    bonds = get_bond_array(sphere_delaunay_bonds_output)
    bond_set = get_bond_set(bonds)
    triangle_list = get_triangle_list(bond_set)

    return np.asarray(triangle_list), np.asarray(list(bond_set))


def write_ply_file(xyz_file_name, xyzs, tris, path):
    file_path = path
    if not file_path.endswith("/"):
        file_path += "/"
    word_split = xyz_file_name.split("_")
    ply_dir = file_path + word_split[0] + "_" + word_split[1]
    if not os.path.exists(ply_dir):
        os.makedirs(ply_dir)
    ply_file_path = ply_dir + "/" + xyz_file_name[:-3] + "ply"
    with open(ply_file_path, "w") as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("comment Created by Ali Farnudi, source file:\n")
        f.write("element vertex {}\n".format(xyzs.shape[0]))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("element face " + str(tris.shape[0]) + "\n")
        f.write("property list uchar uint vertex_indices\n")
        f.write("end_header\n")
        for coords in xyzs:
            f.write("{:.7f} {:.7f} {:.7f}\n".format(coords[0], coords[1], coords[2]))
        for verts in tris:
            f.write("3 {} {} {}\n".format(verts[0], verts[1], verts[2]))


def gen_ply_file(xyz_file_name, args):
    run_command_quiet_verbose(args, "./sphere_delaunay_bonds " + xyz_file_name)

    bond_file_name = xyz_file_name + "l"
    xyzs = np.loadtxt(xyz_file_name)

    tris, bonds = build_triangle_from_bonds(bond_file_name)

    if not args.quiet:
        print_properties(tris, bonds, xyzs)

    tris += -1

    write_ply_file(xyz_file_name, xyzs, tris, args.path)


def gen_bond_info(xyz_file_name, args):
    exe_name = "sphere_delaunay_bonds"
    if not os.path.exists(exe_name):
        run_command(args, "bash sphere_delaunay.sh")
    run_command_quiet_verbose(args, f"./{exe_name} {xyz_file_name}")


def clean_files_with_extention(extention):
    for file in os.listdir():
        if file.endswith(extention):
            os.system(f"rm {file}")


def run_command(args, command):

    if args.verbose:
        os.system(command)
    elif args.quiet:
        os.system(command + " >&- 2>&-")
    else:
        os.system(command + " >&-")


def run_command_quiet_verbose(args, command):

    if args.verbose:
        os.system(command)
    else:
        os.system(command + " >&- 2>&-")


def get_arguments_from_parser():
    parser = argparse.ArgumentParser(
        "ply file generator", description="Generate a random triangulated sphere."
    )
    parser.add_argument(
        "-n",
        help="The number of nodes on the sphere",
        type=int,
        default=1002,
        required=True,
        metavar="",
    )
    parser.add_argument(
        "--max_seed",
        help="maximum number used for the seed. The programmes begins from min_seed and goes up to (not including) this number.",
        type=int,
        required=True,
        metavar="",
    )
    parser.add_argument(
        "--min_seed",
        help="minimum number used for the seed. The programmes begins with this seed number and goes up to (not including) max_seed.",
        type=int,
        default=0,
        metavar="",
    )
    parser.add_argument(
        "--overwrite",
        help="if you want to overwrite xyz files if they exist.",
        action="store_true",
    )
    parser.add_argument(
        "--path",
        help="Specify the path you wish to save the generated ply",
        default="",
    )
    parser.add_argument(
        "-p","--platformID",
        help="VCM platform to run",
        type=int,
        default=2,
    )
    parser.add_argument(
        "-d","--platformDeviceID",
        help="Device ID to use on the platform.",
        type=int,
        default=0,
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-v", "--verbose", action="store_true", help="Display all outputs"
    )
    group.add_argument("-q", "--quiet", action="store_true", help="Silence all outputs")

    return parser.parse_args()


def print_seed(s):
    if args.verbose:
        print("=============================================")
        print("==================== {s} ====================")
        print("=============================================")
    elif args.quiet:
        print(f"S {s}")
    else:
        print(f"----> S {s}")


if __name__ == "__main__":
    args = get_arguments_from_parser()

    for s in range(args.min_seed, args.max_seed):
        print_seed(s)
        xyz_file_name = f"USphere_{args.n}d_s{s}.xyz"
        gen_xyz_command = f"./NodGenClang9 {args.n} {s} {args.platformID} {args.platformDeviceID}"

        if os.path.exists(xyz_file_name):
            if args.overwrite:
                run_command(args, gen_xyz_command)
        else:
            run_command(args, gen_xyz_command)

        gen_bond_info(xyz_file_name, args)
        gen_ply_file(xyz_file_name, args)
    clean_files_with_extention(".eps")
    clean_files_with_extention(".xyzl")
