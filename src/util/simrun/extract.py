import os
import sys
import xml.etree.ElementTree as ET

experiment = sys.argv[1]

print "PAM,RECALL,PRECISION,TIME,MEMORY"
for i in xrange(0, 9):
	id = max(1, i * 3)
	path = os.path.join("result", experiment, str(id), "out.xml")
	tree = ET.parse(path)
	root = tree.getroot()
	now_result = []
	for test in root.findall("homologyTests"):
		now_result.append(test.find("aggregateResults").find("all").attrib["average"])
	rt = [line for line in open(os.path.join("result", experiment, "rt" + str(id) + ".txt"))]
	running_time = float(rt[0].split()[0])
	memory = float(rt[1].split()[0]) / 1000
	recall, precision = now_result
	print ",".join((str(x) for x in (id, recall, precision, running_time, memory)))

