import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import xml.etree.ElementTree as ET

experiment = sys.argv[1]
mypath = join('./result/sibelia/', experiment)
k_value = []
for f in listdir(mypath):
        p = join(mypath, f)
	if isdir(p) and int(f) < 45:
		k_value.append(int(f))

k_value.sort()
for now_k_value in k_value:
	f = str(now_k_value)
	p = join(mypath, str(now_k_value))
        if isdir(p):
		print 'k =', f
		result = dict()
		best_recall = 0
		best_branch = 0
		for f2 in listdir(p):
			p2 = join(p, f2)
			branch = int(f2)
			if branch not in result:
				result[branch] = dict()
			path = os.path.join(p2, "out.xml")
			tree = ET.parse(path)
			root = tree.getroot()
			now_result = []
			for test in root.findall("homologyTests"):
				now_result.append(float(test.find("aggregateResults").find("all").attrib["average"]))
			recall, precision = now_result
			result[branch] = (recall, precision)
			if recall > best_recall:
				best_recall, best_branch = recall, branch
		print "Best recall =", best_recall
		print "Branch =", best_branch
		print "branch,result"
		keys = result.keys()
		keys.sort()
		for key in keys:
			print str(key) + "," + str(result[key][0]) + "/" + str(result[key][1])
		print ""

