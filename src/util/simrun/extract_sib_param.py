import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import xml.etree.ElementTree as ET


print "Experiment,Best k,Best branch,Recall,Precision,Sibelia RT,Sibelia memory,spoa RT,spoa memory"

for experiment_n in xrange(0, 9):
        experiment = str(max(experiment_n * 3, 1))
	mypath = join('./result/sibelia/', experiment)
	k_value = []
	for f in listdir(mypath):
        	p = join(mypath, f)
		if isdir(p) and int(f) < 45:
			k_value.append(int(f))
	k_value.sort()
	param_file = open("sib_param" + experiment + ".txt", "w")
	best_k,best_global_branch,best_global_recall,best_global_precision=(0,0,0,0)
	for now_k_value in k_value:
		f = str(now_k_value)
		p = join(mypath, str(now_k_value))
        	if isdir(p):
			print >> param_file, 'k =', f
			result = dict()
			best_recall = 0
			best_branch = 0
			for f2 in listdir(p):
				p2 = join(p, f2)
				branch = int(f2)
				path = os.path.join(p2, "out.xml")
				tree = ET.parse(path)
				root = tree.getroot()
				now_result = []
				for test in root.findall("homologyTests"):
                                        val = float(test.find("aggregateResults").find("all").attrib["average"]) * 100
                                        now_result.append(int(round(val)))
				recall, precision = now_result
				result[branch] = (recall, precision)
				rt_memory = []
				for line in open(os.path.join(p2,  "rt" + str(experiment) + ".txt")):
					line = line.strip().split()
					rt_memory.append(int(round(float(line[1]))))
					rt_memory.append(int(line[3]) / 1000)
				if recall >= best_recall:
					best_recall, best_branch = recall, branch
				if recall >= best_global_recall:
					best_k,best_global_branch,best_global_recall,best_global_precision = now_k_value, branch,  recall, precision
					best_rt_memory = rt_memory
			print >> param_file, "Best recall =", best_recall
			print >> param_file, "Branch = ", best_branch
			print >> param_file, "Branch,Recall/Precision"
			keys = result.keys()
			keys.sort()
			for key in keys:
				print >> param_file, str(key) + "," + str(result[key][0]) + "/" + str(result[key][1])
			print >> param_file, ""
	print ','.join(str(x) for x in [experiment,best_k,best_global_branch,best_global_recall,best_global_precision] + best_rt_memory)


