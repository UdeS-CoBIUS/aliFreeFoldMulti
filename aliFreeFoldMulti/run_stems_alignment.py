#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Libraries and files
from align_seq_struct3 import main

# Constant variables
OUTPUTS = "../aliFreeFoldMulti_Results/"


def run_stems_alignment(method_name, output_folder_path, is_verbose):
    """ Run the stems alignment using the stems alignment method

    :param method_name: Method name for the stems alignment
    :type method_name: str
    :param output_folder_path: Folder path for the output files
    :type output_folder_path: str
    :param is_verbose: Display the outputs of stems alignment
    :type is_verbose: bool

    :return: None

    """
    rep_file = output_folder_path + "rep.db"
    subopt_file = output_folder_path + "subopt.db"

    if is_verbose:
        stems_alignment_outputs = main(rep_file, subopt_file, withWobble=True,
                                       method=method_name, oldRun=False,
                                       extendRef=True, extendRNA=True,
                                       verbose=True)
    else:
        stems_alignment_outputs = main(rep_file, subopt_file, withWobble=True,
                                       method=method_name, oldRun=False,
                                       extendRef=True, extendRNA=True,
                                       verbose=False)

    stems_alignment_outputs_file = "{}{}-alignment.db".format(
        output_folder_path, method_name)

    stems_alignment_file = open(stems_alignment_outputs_file, "w")

    for output_tuple in stems_alignment_outputs:
        stems_alignment_file.write(">{}\n".format(output_tuple[0]))
        stems_alignment_file.write("{}\n".format(output_tuple[1]))
        stems_alignment_file.write("{}\n".format(output_tuple[2]))

    stems_alignment_file.close()
