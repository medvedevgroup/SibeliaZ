import csv

param_table = open("param_table.txt", "w")
rt_memory_table = open("rt_memory_table.txt", "w")

def mem(val):
	return " (" + val + ")"

with open('cactus.csv') as cactus_csvfile, open('sibelia.csv') as sibelia_csvfile, open('lsibelia.csv') as lsibelia_csvfile:
	cactus = csv.DictReader(cactus_csvfile)
        sibelia = csv.DictReader(sibelia_csvfile)
        lsibelia = csv.DictReader(lsibelia_csvfile)
	for cactus_row in cactus:
		sibelia_row = sibelia.next()
		lsibelia_row = lsibelia.next()
		experiment = sibelia_row["Experiment"]
		assert cactus_row["Experiment"] == experiment
		assert lsibelia_row["Experiment"] == experiment
		lsibelia_data = [lsibelia_row["TwoPaCo RT"] +  mem(lsibelia_row["TwoPaCo memory"]), lsibelia_row["L-Sibelia RT"] +  mem(lsibelia_row["L-Sibelia memory"]), lsibelia_row["spoa RT"] +  mem( lsibelia_row["spoa memory"])]
		sibelia_data = [sibelia_row["Sibelia RT"] +  mem(sibelia_row["Sibelia memory"]), sibelia_row["spoa RT"] +  mem(sibelia_row["spoa memory"])]
		cactus_data = [cactus_row["Cactus RT"] +  mem(cactus_row["Cactus memory"])]
		print >> rt_memory_table, "&".join([experiment] + lsibelia_data + sibelia_data + cactus_data) + "\\\\"
		lsibelia_param = [lsibelia_row["Best k"], lsibelia_row["Best branch"], lsibelia_row["Best flank"]]
		sibelia_param = [sibelia_row["Best k"], lsibelia_row["Best branch"]]
		print >> param_table, "&".join([experiment] + lsibelia_param + sibelia_param) + "\\\\"


