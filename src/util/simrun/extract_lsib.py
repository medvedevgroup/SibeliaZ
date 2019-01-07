import os
import sys
import xml.etree.ElementTree as ET

experiment = sys.argv[1]

print "PAM,RECALL,PRECISION,TIME1,MEMORY1,TIME2,MEMORY2,TIME3,MEMORY3"
for i in xrange(0, 9):
	id = max(1, i * 3)
	path = os.path.join("result", experiment, str(id), "out.xml")
	tree = ET.parse(path)
	root = tree.getroot()
	now_result = []
	for test in root.findall("homologyTests"):
		now_result.append(test.find("aggregateResults").find("all").attrib["average"])
	rt = [line for line in open(os.path.join("result", experiment, "rt" + str(id) + ".txt"))]
	recall, precision = now_result
	output = [id, recall, precision]
	for j in xrange(0, 3):
		running_time = float(rt[j * 2].split()[0])
		memory = float(rt[j * 2 + 1].split()[0]) / 1000
		output.append(running_time)
		output.append(memory)
	print ",".join((str(x) for x in output))

