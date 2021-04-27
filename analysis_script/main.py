"""
Author: David Ho
Institute: National Tsing Hua university, Department of Physics, Hsinchu, Taiwan 
Mail: davidho@gapp.nthu.edu.tw
"""
#Import packages
from script import cutflow, parse, chi2, chi2_from_npz, fitting 
import h5py, sys, traceback, os, tqdm
from argparse import ArgumentParser
import multiprocessing as mp

def main():
    
    parser = ArgumentParser()
    parser.add_argument("-u", "--usage", dest="usage", help="Define the purpose to run.")
    parser.add_argument("-i", "--input-file", dest="input",  help="Input file directory.")
    parser.add_argument("-o", "--output-file", dest="output",  help="Output file directory.")
    parser.add_argument("-m", "--model", dest="model", help="Select a model to parse. Usable mode: ttbar, ttH, four_top, bkg")
    parser.add_argument("-c", "--config", dest="config", help="The directory configuration file for cutflow.")
    parser.add_argument("-s", "--single", dest="single", default=True, help="Determin is dealing with single file or not. If not dealing with single file, please input the directory of root files.")
    parser.add_argument("-p", "--num_process", dest="process", default=1, type=int, help="Number of extra process for accelerating speed.")
    parser.add_argument("-e", "--exrta_option", dest="extra", default="normal", help="Extra option for used.")
    parser.add_argument("-t", "--target_file", dest="target", help="target file for computing purity")
    parser.add_argument("-g", "--shower-generator", dest="generator", default="py8", help="type of generator.")

    args = parser.parse_args()
    
    if args.single == str(True):
        args.single = True
    else: 
        args.single = False

    if args.usage == "cutflow":
        cutflow(args.input, args.output, args.model, args.config, args.single)
    elif args.usage == "parse":
        parse(args.input, args.output, args.model, args.process, args.generator, args.single)
    elif args.usage == "chi2":
        chi2(args.input, args.output, args.model, args.process, args.extra, args.generator, args.single)
    elif args.usage == "chi2_from_npz":
        chi2_from_npz(args.input, args.output, args.model, args.process, args.extra, args.single)
    # elif args.usage == "purity":
    #     purity(args.input, args.target, args.output, args.model, args.single)
    elif args.usage == "fitting":
        fitting(args.input, args.output, args.model, args.single)
    # elif args.usage == 'background':
    #     background(args.input, args.output, args.model, args.single, args.process, args.extra, args.generator)
    else: 
        print("Please select a correct usage.\n1. cutflow\n2. parse")


if __name__ == '__main__':
    main()
    
