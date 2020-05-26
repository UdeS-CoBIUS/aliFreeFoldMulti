#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Name: aliFreeFoldMulti program
Author: Marc-Andre Bossanyi
Description: Alignment-free approach of secondary structure prediction
"""

# Libraries and files
from time import time

from parse_arguments import parse_arguments
from postprocessing import remove_files, replace_files_name
from run_alifreefold import run_alifreefold_using_centroid, \
    run_alifreefold_using_subopts
from run_stems_alignment import run_stems_alignment

# Constant variables
ALL_METHODS = ["first", "last", "best"]


def main():
    """ Main function for aliFreeFoldMulti

    :return: None

    """
    input_file_name, method_name, output_folder_path, is_centroids, \
        is_subopt, is_verbose = parse_arguments()

    initial_time = time()

    # == ALIFREEFOLD FOR CENTROID ==
    print("Running aliFreeFoldMulti ...")
    run_alifreefold_using_centroid(input_file_name, output_folder_path,
                                   is_verbose)

    # == STEMS ALIGNMENT ==
    if method_name == "all":
        for method in ALL_METHODS:
            run_stems_alignment(method, output_folder_path, is_verbose)
    else:
        run_stems_alignment(method_name, output_folder_path, is_verbose)

    # == ALIFREEFOLD FOR SUBOPTS ==
    if method_name == "all":
        for method in ALL_METHODS:
            run_alifreefold_using_subopts(method, output_folder_path,
                                          is_verbose)
    else:
        run_alifreefold_using_subopts(method_name, output_folder_path,
                                      is_verbose)

    # == COMPUTE TIME ==
    total_time = time() - initial_time

    # == POSTPROCESSING FILES ==
    replace_files_name(output_folder_path)
    remove_files(output_folder_path, is_centroids, is_subopt)

    # == END ==
    print("Total time: {}s".format(total_time))


if __name__ == "__main__":
    main()
