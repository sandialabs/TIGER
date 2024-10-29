import networkx as nx
import csv
import os
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
import csv
import configargparse
import sys

def network(ANI):
	'''
	given an ANI file, output:
	[name1, name2, ANIdist, AF]
	'''
	with open(ANI) as g:
		#print 'played it till my fingers bled'
		reader = csv.reader(g, delimiter="\t")
		y = list(reader)
		toggle='no'
		#print 'i got my first real six string'
		output=[]
		count=0
		for line in y:
			count+=1
			#print line
			output.append([line[0],line[1],1-(float(float(line[2])/100)),float(line[3])/float(line[4])])
			#output.append([line[0],line[1],1-(float(float(line[2])/100))])
			#print count
	return output

def getAve(islandComp,dict):
    '''
    given a line from the ANI list, get the average scores
    don't worry about done list, multiprocessing should speed this up enough for it not to matter
    '''
    if islandComp[0] in dict[islandComp[1]]:
        return [islandComp[0],islandComp[1],(float(islandComp[2])+float(dict[islandComp[1]][islandComp[0]][0]))/float(2),(float(islandComp[3])+float(dict[islandComp[1]][islandComp[0]][1]))/float(2)]
    else:
        return [islandComp[0],islandComp[1],float(islandComp[2])+1/float(2),float(islandComp[3])+0/float(2)]



if __name__ == '__main__':

	if len (sys.argv) < 5 :
		print("Usage: python fastANIclust.py -i input_edge_file -o output_prefix -d dashing? -an ANI_cutoff -af AF_cutoff\nNumber of arguements is " + str(len(sys.argv)))
		sys.exit (1)
	p = configargparse.ArgParser(description='this script makes network files from distance matrices')
	p.add('-i', required=True, help='input file',dest='file')
	p.add('-o', required=True, help='output prefix',dest='out') 
	p.add('-d', required=False, help='dashing',dest='dash', action='store_false') 
	p.add('-an', required=True, help='ANI cutoff (higher is more stringent)',dest='ani') 
	p.add('-af', required=True, help='AF cutoff (higher is more stringent)',dest='af')

	#print 'we have started'
	args=p.parse_args()

    #file="Islands_Out70_FastANI_all.out"
    #outputfile='Islands_Out70_FastANI_all_network.out'

    
    final=network(args.file)
    
    nameOfFixedFile = str(args.out)+'.proc'
    
    with open(args.out,"w") as g:
        wr = csv.writer(g)
        wr.writerows(final)

    with open(nameOfFixedFile,"r") as f:
        reader = csv.reader(f, delimiter=",")
        x = list(reader)

    #version 3 with dictionaries
    ANIdict={}
    for entry in x:
        if entry[0] in ANIdict:
            ANIdict[entry[0]][entry[1]]=[entry[2],entry[3]]
        else:
            ANIdict[entry[0]]={}
            ANIdict[entry[0]][entry[1]]=[entry[2],entry[3]]


    end=[]
    test=0
    done={}
    test2=0
    test3=0
    test4=0
    test5=0
    for entry in x:
        test+=1
        if test%10000==0:
            print str(datetime.now())
        if entry[1] in done:
            if entry[0] in done[entry[1]]:
                #print 'hello'
                test2+=1
                continue
            else:
                end.append(getAve(entry,ANIdict))
                test3+=1
        elif entry[0] in done:
            end.append(getAve(entry,ANIdict))
            done[entry[0]].append(entry[1])
            test4+=1
        else:
            end.append(getAve(entry,ANIdict))
            done[entry[0]]=[entry[1]]
            test5+=1

    outfile=str(args.out)+'_ave'

    with open(outfile,"w") as f:
        wr = csv.writer(f,delimiter=",")
        wr.writerows(end)
            


    ##remove AF value, and filter, from file
    with open(outfile,"r") as f:
        reader = csv.reader(f, delimiter=",")
        x = list(reader)
        end=[]
        for line2 in x:
            #print(line2)
            if float(line2[2]) < 0.05 and float(line2[3]) > 0.90:
                if float(ANIdict[line2[0]][line2[1]][0]) < 0.05 and float(ANIdict[line2[1]][line2[0]][0]) < 0.05 and float(ANIdict[line2[0]][line2[1]][1]) > 0.9 and float(ANIdict[line2[1]][line2[0]][1]) > 0.9:
                #print 'yes'
                    if line2[0]<line2[1]:
                        end.append([line2[0],line2[1],line2[2]])
                    else:
                        end.append([line2[1],line2[0],line2[2]])

    out2=str(outfile)+'_ANI'+str(args.ani)+'_AF'+str(args.af)

    with open(out2,"w") as f:
        wr = csv.writer(f,delimiter=",")
        wr.writerows(end)

    ##load for graphing and clustering
    print('Files all pre-processed. Floyd Warshall launching.')
    G=nx.read_edgelist(out2,delimiter=',', nodetype=str,
      data=(('weight',float),))

    #can we get connected components
    #test=sorted(nx.connected_components(G), key = len, reverse=True)
    #set(test[0]).intersection(test[1])

    rk=nx.floyd_warshall_numpy(G, nodelist=G.nodes(), weight='weight')
    rk[rk == inf]=1
    print('Floyd Warshall complete. DBSCAN launching.')
    rk2=DBSCAN(eps=1, min_samples=2).fit_predict(rk)
    print('DBSCAN complete. Writing files.')
    clusterout={}
    for a in range(0,len(rk2)):
        if rk2[a] in clusterout:
            clusterout[rk2[a]].append(G.nodes()[a])
        else:
            clusterout[rk2[a]]=[G.nodes()[a]]

    ##remove path from clusterout
    #for k,v in clusterout.iteritems():
    #	newlist=[]
    #	for a in v:
    #		newa=a.replace('/home/rkrishn/projects/islands/Raga/FastANI/out70/','')
    #		newlist.append(newa)
    #	clusterout[k]=newlist

    clusterName=str(out2)+'_clusters.txt'
    clusterPickle=str(out2)+'_clusters.pickle'

    ##save clusters
    with open(clusterName, 'wb') as handle:
        pickle.dump(clusterout, handle, protocol=pickle.HIGHEST_PROTOCOL)

    ##save clusters as txt file
    with open(clusterPickle, 'w') as f:
        for k,v in clusterout.iteritems():
            if k == -1:
                pass
            else:
                f.write('\n%s\t' % (k))
                for a in v:
                    f.write('%s\t' % (a))
                    
    print('Files written')