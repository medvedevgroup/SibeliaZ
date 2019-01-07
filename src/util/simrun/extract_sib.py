import os
import sys
import xml.etree.ElementTree as ET

experiment = sys.argv[1]

print "PAM,RECALL,PRECISION"
for i in xrange(0, 9):
	id = max(1, i * 3)
	path = os.path.join("result", experiment, str(id), "out.xml")
	tree = ET.parse(path)
	root = tree.getroot()
	now_result = []
	for test in root.findall("homologyTests"):
		now_result.append(test.find("aggregateResults").find("all").attrib["average"])
	recall, precision = now_result
	print ",".join((str(x) for x in (id, recall, precision)))

