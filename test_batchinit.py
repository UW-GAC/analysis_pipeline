#! /usr/bin/env python3
import os
import sys

from       argparse import ArgumentParser

try:
    import boto3
    import batchInit
except ImportError:
    print ("AWS batch not supported.")



# parse input
parser = ArgumentParser( description = "script to test batchinit" )
parser.add_argument( "-n", "--name",
                     help = "Name of the compute environment to create [required]" )
parser.add_argument( "-p", "--profile", default = "uw",
                     help = "aws profile")
parser.add_argument( "-q", "--queue",
                     help = "batch queue name [required]")
parser.add_argument( "-i", "--instancetypes", default = "m5.2xlarge, m5.4xlarge, c5.4xlarge",
                     help = "batch queue name [required]")
parser.add_argument( "--pricing", default = "EC2",
                     help = "instance pricing  [EC2 | SPOT]")
parser.add_argument( "-t", "--tag", default = "111222333",
                     help = "aws tag")
parser.add_argument( "-D", "--Debug", action="store_true", default = False,
                     help = "Turn on verbose output [default: False]" )
args = parser.parse_args()
name = args.name
profile = args.profile
queue = args.queue
instancetypes = args.instancetypes
pricing = args.pricing
tag = args.tag
debug = args.Debug

if name == None:
    print("Name (-n) is required.")
    sys.exit(2)

if queue == None:
    print("queue (-q) is required.")
    sys.exit(2)

try:
    session = boto3.Session(profile_name = profile)
    batchC = session.client('batch')

except Exception as e:
    print('boto3 session or client exception ' + str(e))
    sys.exit(2)

batchInit.createEnv(batchC, queue, pricing, instancetypes, name, tag, profile, verbose=debug)
