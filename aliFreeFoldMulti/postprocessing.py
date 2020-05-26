#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Librairies and files
from os import listdir, remove, replace
from os.path import isfile, join


def replace_files_name(output_folder_path):
    """ Replace the name of the files from the output folder

    :param output_folder_path: Folder path for the output files
    :type output_folder_path: str

    :return: None

    """
    # Replace the name of the file "strat1.db"
    old_file_path = "{}strat1.db".format(output_folder_path)
    new_file_path = "{}centroids.db".format(output_folder_path)
    replace(old_file_path, new_file_path)

    # Replace the name of the files "-alignment.db"
    files = [f for f in listdir(output_folder_path)
             if isfile(join(output_folder_path, f))
             if f.strip().split("-")[-1] == "alignment.db"]

    for file in files:
        old_file_path = "{}{}".format(output_folder_path, file)
        new_file_path = "{}{}-stems_alignment.db".format(output_folder_path,
                                                         file.split("-")[0])
        replace(old_file_path, new_file_path)

    # Replace the name of the files "-strat_a.db"
    files = [f for f in listdir(output_folder_path)
             if isfile(join(output_folder_path, f))
             if f.strip().split("-")[-1] == "strat_a.db"]

    for file in files:
        old_file_path = "{}{}".format(output_folder_path, file)
        new_file_path = "{}{}-closest_subopt.db".format(output_folder_path,
                                                        file.split("-")[0])
        replace(old_file_path, new_file_path)


def remove_files(output_folder_path, is_centroids, is_subopt):
    """ Delete the files not useful for the user

    :param output_folder_path: Folder path for the output files
    :type output_folder_path: str
    :param is_centroids: Centroids option is selected
    :type is_centroids: bool
    :param is_subopt: Closest subopt option is selected
    :type is_subopt: bool

    :return: None

    """
    # Delete the uploaded file ".fa", if it exists
    file = [f for f in listdir(output_folder_path)
            if isfile(join(output_folder_path, f))
            if f.strip().split(".")[-1] == "fa"
            or f.strip().split(".")[-1] == "fasta"]

    if len(file) > 0:
        file_to_remove = "{}{}".format(output_folder_path, file)
        remove(file_to_remove)

    # Delete the files with the format ".csv"
    files = [f for f in listdir(output_folder_path)
             if isfile(join(output_folder_path, f))
             if f.strip().split(".")[-1] == "csv"]

    for file in files:
        file_to_remove = "{}{}".format(output_folder_path, file)
        remove(file_to_remove)

    # Delete the file "rep.db"
    file_to_remove = "{}rep.db".format(output_folder_path)
    remove(file_to_remove)

    # Delete the file "strat2.db"
    file_to_remove = "{}strat2.db".format(output_folder_path)
    remove(file_to_remove)

    # Delete the file "subopt.db"
    file_to_remove = "{}subopt.db".format(output_folder_path)
    remove(file_to_remove)

    # If the option "centroids" is not selected, delete
    if not is_centroids:
        file_to_remove = "{}centroids.db".format(output_folder_path)
        remove(file_to_remove)

    # If the option "subopt" is not selected, delete
    if not is_subopt:
        files = [f for f in listdir(output_folder_path)
                 if isfile(join(output_folder_path, f))
                 if f.find("_") != -1]

        for file in files:
            file_to_remove = "{}{}".format(output_folder_path, file)
            remove(file_to_remove)

