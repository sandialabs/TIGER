#wrapper
import subprocess
import os
import networkx as nx
import csv
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.pylab as plt2
import numpy as np
import scipy
import pickle
import sklearn
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import pandas
from collections import Counter
from numpy import inf
from os import walk
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs
import json
from datetime import datetime
import configargparse


cwd = os.getcwd()
	
if len (sys.argv) < 11 :
	print "Usage: python islandCluster.py -i input_edge_file -an ANI_cutoff -af AF_cutof -p path_for_results (default is current path) -x prefix_for_output\nNumber of arguements is " + str(len(sys.argv))
	sys.exit (1)
print 'very beginning'
p = configargparse.ArgParser(description='this script makes network files from distance matrices')
p.add('-i', required=True, help='input file (multifasta of islands) - use full path',dest='file') 
p.add('-an', required=True, help='ANI cutoff (higher is more stringent)',dest='ani') 
p.add('-af', required=True, help='AF cutoff (higher is more stringent)',dest='af')
p.add('-p', required=True, help='path for results', dest='path', default=cwd)
p.add('-x', required=True, help='prefix for results', dest='prefix')
p.add('-n', required=True, help='number of processes to spawn', dest='num', default=1)
args=p.parse_args()

##Jan 2019 FastANI
#make folder
folder=str(args.path)+'/'+str(args.prefix)+'/'
subprocess.Popen("mkdir " + folder, shell=True)
#move fasta
nameOfFile=args.file
nameLite=nameOfFile.split('/')[-1]
nameMove=str(folder)+str(nameLite)
subprocess.Popen("scp " +  nameOfFile + ' ' +  nameMove, shell=True)
#copy file and remove anything after space in header
subprocess.Popen("sed -i '/^>/ s/ .*//' " + nameMove, shell=True)
subprocess.Popen("sed -i 's,/,_,g' " + nameMove, shell=True)

#split into multi fasta

subprocess.Popen("cat " + nameMove + " | awk '{if (substr($0, 1, 1)=='>') {filename=(substr($0,2) '.fa')} print $0 > filename}'", shell=True)
subprocess.Popen("mv " + nameMove+ " ../", shell=True)

##make list of files
locName=str(nameMove)+'_locations.txt'
subprocess.Popen("ls -d -1 $PWD/** > " + locName, shell=True)

numproc=args.num
outANI=str(args.path)+str(nameLite)+'_FastANI.out'
subprocess.Popen("fastANI -k 15 -t " + numproc + " --fragLen 500 --minFrag 2 --ql " + locName + " --rl " + locName + " -o " + outANI, shell=True)

##remove path from file
subprocess.Popen("sed -i 's," + folder + ",,g' " + outANI, shell=True)

##run clustering
subprocess.Popen("python fastANIclust.py -i input_edge_file -o output_file -an " + args.ani + " -af " + args.af, shell=True)

