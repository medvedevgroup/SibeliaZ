import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import xml.etree.ElementTree as ET

print "Experiment,Recall,Precision,Cactus RT,Cactus memory,Conversion RT,Conversion memory"

for experiment_n in xrange(0, 9):
        experiment = str(max(experiment_n * 3, 1))
	mypath = join('./result/cactus/', experiment)
	path = os.path.join(mypath, "out.xml")
	tree = ET.parse(path)
	root = tree.getroot()
	now_result = []
	for test in root.findall("homologyTests"):
		val = float(test.find("aggregateResults").find("all").attrib["average"]) * 100
		now_result.append(int(round(val)))
	rt_memory = []
	recall, precision = now_result
	for line in open(os.path.join(mypath, "rt" + str(experiment) + ".txt")):
		line = line.strip().split()
		rt_memory.append(int(round(float(line[1]))))
		rt_memory.append(int(line[3]) / 1000)
	print ','.join(str(x) for x in [experiment,recall,precision] + rt_memory)


