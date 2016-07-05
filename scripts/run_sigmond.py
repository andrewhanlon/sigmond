#!/usr/bin/env python

import argparse
import xml.etree.ElementTree as xml_handler
import subprocess
import os
import sigmondInput
import sys

def main():

  # Get command line arguments

  parser = argparse.ArgumentParser(description="Run Sigmond Analysis")

  parser.add_argument("-i", "--interactive", action="store_true", required=False,
                      help="Start interactive mode (Not currently implemented)")
  parser.add_argument("--do-checks", action="store_true", required=False,
                      help="Do the DoChecks task")
  parser.add_argument("--merge", action="store_true", required=False,
                      help="Merge all channels into one XML input file")
  parser.add_argument("--exclude-dirs", nargs='+', required=False, default=[],
                      help="Specify a list of directories to exclude")
  parser.add_argument("-d", "--dir", type=str, required=False, default=".",
                      help="Specify base directory. Defaults to current directory")
  parser.add_argument("--list-subdirs", action="store_true", required=False,
                      help="List all subdirectories that will be searched")
  parser.add_argument("--list-ops", action="store_true", required=False,
                      help="List all operators found")
  parser.add_argument("--output-ops", type=str, required=False, metavar="FILE NAME",
                      help="Output the list of operators to a text file that")
  parser.add_argument("--specify-ops", type=str, required=False, metavar="FILE NAME",
                      help="Specify which operators you would like to include. Pass a text file with each operator on a separate line")
  parser.add_argument("--herm-corr", action="store_true", required=False, default=False,
                      help="Set specifications as a Hermitian correlator matrix")
  parser.add_argument("-x", "--execute", action="store_true", required=False, default=False,
                      help="Execute sigmond after generating the input XML. The default is to create the XML only")
  parser.add_argument("-n", "--name", type=str, required=False, default="",
                      help="Specify name of Project. Default is nothing. Also used to name input XML and logfile")
  parser.add_argument("-e", "--echo", action="store_true", required=False,
                      help="Specify whether the input XML should be echoed in the log file")
  parser.add_argument("-v", "--verbose", action="store_true", required=False,
                      help="Add verbosity tags")
  parser.add_argument("--printXML", type=str, required=False,
                      choices=['MCValues','MCBootstraps','MCJackknives','MCHistogram','MCBootstrapHistogram','MCJackknifeHistogram'],
                      help="Performs the PrintXML task for every MCObservable")
  parser.add_argument("--print-corr", action="store_true", required=False,
                      help="Performs the PrintXML task of type TemporalCorrelator for every MCObservable")
  parser.add_argument("--print-energy", action="store_true", required=False,
                      help="Performs the PrintXML task of type EffectiveEnergy for every MCObservable")
  parser.add_argument("--real", action="store_true", required=False,
                      help="Includes <Arg>Re</Arg> tags")
  parser.add_argument("--imag", action="store_true", required=False,
                      help="Includes <Arg>Im</Arg> tags")
  parser.add_argument("--bins", type=int, required=False, default=40,
                      help="Specifiy number of bins")
  parser.add_argument("--do-fits", action="store_true", required=False,
                      help="Performs fits")
  parser.add_argument("--do-fits-plots", action="store_true", required=False,
                      help="Performs fits and plots them")
  parser.add_argument("--bootstrap", action="store_true", required=False, default=False,
                      help="Use bootstrapping. Default is jackknife")
  parser.add_argument("-r", "--read", action="store_true", required=False,
                      help="Read and analyze the results")
  parser.add_argument("--laph-query", type=str, required=False, metavar="LAPH_QUERY", default="laph_query",
                      help="Specify LapH query executable")
  parser.add_argument("--sigmond", type=str, required=False, metavar="SIGMOND", default="sigmond",
                      help="Specify Sigmond executable")

  args = parser.parse_args()

  # No spaces in the name
  name = args.name.replace(' ', '')
  name = name.replace('\t', '')
  name = name.replace('\n', '')
  name = name.replace('\r', '')

  sig_inputs = []

  # Find data directories
  for root, dirs, files in os.walk(args.dir):
    
    if not set(root.split('/')).isdisjoint(args.exclude_dirs):
      continue

    if args.list_subdirs:
      print(root)
      continue

    if contains_laph_data(args.laph_query, root):
      print("Adding data from: ", root)
      sig_inputs.append(sigmondInput.SigmondInput(name, root, args.echo, args.bootstrap, args.laph_query))

  if args.merge:
    merged_sig = sum(sig_inputs)
    sig_inputs[:] = []
    sig_inputs.append(merged_sig)

  for sig in sig_inputs:
    
    if args.herm_corr:
      sig.set_herm()
    
    if (args.output_ops):
      sig.print_ops(sig.channel_name + "_" + args.output_ops)
  
    if args.list_ops:
      sig.print_ops()

    # Tasks

    if args.printXML:
      if args.real:
        sig.printXML(args.printXML, True, args.bins, args.verbose)

      if args.imag:
        sig.printXML(args.printXML, False, args.bins, args.verbose)

      if not args.imag and not args.real:
        sig.printXML(args.printXML, True, args.bins, args.verbose)
        sig.printXML(args.printXML, False, args.bins, args.verbose)


    if args.print_corr:
      if args.real:
        sig.printCorr(True)

      if args.imag:
        sig.printCorr(False)

      if not args.imag and not args.real:
        sig.printCorr(True)
        sig.printCorr(False)

    if args.print_energy:
      if args.real:
        sig.printEnergy(True)

      if args.imag:
        sig.printEnergy(False)

      if not args.imag and not args.real:
        sig.printEnergy(True)
        sig.printEnergy(False)

    if args.do_fits:
      sig.do_fits()

    if args.do_checks:
      sig.do_checks()

    # Write to file
    sig.write()
  
    if args.execute:
      sig.execute(args.sigmond)

    if args.read:
      sig.read()


def contains_laph_data(laph_query, direc):
  
  if not os.path.isdir(direc):
    print("WARNING: Invalid directory. No LapH data here\n")
    return False

  for f in os.listdir(direc):
    
    full_file = os.path.join(direc, f)

    if not os.path.isfile(full_file):
      continue

    result = subprocess.check_output([laph_query, full_file]).strip().decode()

    if result == "This is a LapH correlator file" or result == "This is a LapH VEV file":
      return True

  return False

if __name__ == "__main__":

  main()
