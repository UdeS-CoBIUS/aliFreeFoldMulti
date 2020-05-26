#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Libraries and files
import argparse

# Constant variables
METHOD_CHOICES = ["first", "last", "best", "all"]
OPTION_CHOICES = ["true", "false"]


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("i", help="input file name")
    parser.add_argument("m", type=str, choices=METHOD_CHOICES,
                        help="input method name")
    parser.add_argument("o", type=str, help="output folder")
    parser.add_argument("c", type=str, choices=OPTION_CHOICES,
                        help="is 'centroids' option selected")
    parser.add_argument("s", type=str, choices=OPTION_CHOICES,
                        help="is 'subopt' option selected")
    parser.add_argument("--verbose", help="increase output verbosity")

    args = parser.parse_args()

    input_file_name = args.i
    method_name = args.m
    output_folder_path = args.o
    is_centroids = args.c
    is_subopt = args.s
    is_verbose = args.verbose

    if is_centroids == "true":
        is_centroids = True
    else:
        is_centroids = False

    if is_subopt == "true":
        is_subopt = True
    else:
        is_subopt = False

    if not is_centroids and not is_subopt:
        print("ERROR: You must select at least one option !")
        print("\nEnd of aliFreeFoldMulti")
        exit()

    return input_file_name, method_name, output_folder_path, is_centroids, \
           is_subopt, is_verbose
