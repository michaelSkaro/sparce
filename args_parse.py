import os
import re
import sys
import argparse
import subprocess



def parse_cla(args):
    '''
    Parses C ommand L ine A rguments.
    Parameters
    ----------
    args: list
        List of command line options followed by arguments.
    
    Returns
    -------
    args: dict
        Parsed arguments in key-value pairs.
    '''
    
    parser = argparse.ArgumentParser(
        description='A supervisory script to run all data through the'
                    + ' feature selection pipeline.'
    )
    required = parser.add_argument_group('required args')
    required.add_argument('-x', '--input', required=True,
                          help='Path to input matrix for features. A pandas dataframe or csv file.')
    required.add_argument('-y', '--target', required=True,
                          help='Target variable.')
    required.add_argument('-f', '--top_n', required=True,
                          help='We will grade each feature for you. This will tailor the filter to select the top n features. This can be ran in series or blocks to recursively eliminate redundant features.')
    required.add_argument('-j', '--number_of_jobs', required=True,
                          help='The number of CPUs to dedicate to the feature selection, default == 0')

    return vars(parser.parse_args(args))
