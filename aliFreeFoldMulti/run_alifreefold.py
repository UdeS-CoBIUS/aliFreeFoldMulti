#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Libraries and files
from os import replace
import platform
from subprocess import call, DEVNULL

# Constant variables
ALIFREEFOLD_PATH = "../alifreefold-master/bin/Release/alifreefold"
ALL_OPTIONS = " -f abcdefghstu"
USING_CENTROID = " -t full"
USING_SUBOPTS = " -t partial"
OUTPUTS = "../aliFreeFoldMulti_Results/"


def run_alifreefold_using_centroid(input_file_name, output_folder_path,
                                   is_verbose):
    """ Run aliFreeFold using the centroid method

    :param input_file_name: File name ".fa"
    :type input_file_name: str
    :param output_folder_path: Folder path for the output files
    :type output_folder_path: str
    :param is_verbose: Display the outputs of aliFreeFold
    :type is_verbose: bool

    :return: None

    """
    if platform.system() == "Windows":
        open_bash = "bash -c cd ./;"
    else:
        open_bash = ""

    command = open_bash + ALIFREEFOLD_PATH + USING_CENTROID \
                        + " -i " + input_file_name \
                        + " -o " + output_folder_path + ALL_OPTIONS

    if is_verbose:
        call(command, shell=True)
    else:
        call(command, shell=True, stdout=DEVNULL)


def run_alifreefold_using_subopts(method_name, output_folder_path,
                                  is_verbose):
    """ Run aliFreeFold using the closest 25-subopt structures method

    :param method_name: Method name for the stems alignment
    :type method_name: str
    :param output_folder_path: Folder path for the output files
    :type output_folder_path: str
    :param is_verbose: Display the outputs of aliFreeFold
    :type is_verbose: bool

    :return: None

    """
    stems_alignment_file_path = output_folder_path + method_name + "-alignment.db"
    subopt_file_path = output_folder_path + "subopt.db"
    weights_file_path = output_folder_path + "weights.csv"
    w_nm_file_path = output_folder_path + "w_nmRepOfSS.csv"

    input_files = ",".join([stems_alignment_file_path, subopt_file_path,
                            weights_file_path, w_nm_file_path])

    if platform.system() == "Windows":
        open_bash = "bash -c cd ./;"
    else:
        open_bash = ""

    command = open_bash + ALIFREEFOLD_PATH + USING_SUBOPTS \
                        + " -i " + input_files + " -o " + output_folder_path

    if is_verbose:
        call(command, shell=True)
    else:
        call(command, shell=True, stdout=DEVNULL)

    old_file_path = "{}strat_a.db".format(output_folder_path)
    new_file_path = "{}{}-strat_a.db".format(output_folder_path, method_name)

    replace(old_file_path, new_file_path)
