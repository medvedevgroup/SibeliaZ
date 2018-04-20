template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "%s_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "strand" />
		<Filter name = "biotype" value = "protein_coding"/>
	</Dataset>
</Query>'''

hmlg_template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "%s_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "%s_homolog_ensembl_gene" />
	</Dataset>
</Query>'''

paralog_template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "%s_gene_ensembl" interface = "default" >
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "%s_paralog_ensembl_gene" />
        </Dataset>
</Query>'''


strain = []
hmlg = []

for line in open("mice.txt"):
	s = line.strip().split()[0]
	strain.append(s)

for s in strain:
	file = open("./gene/" + s + "_gene.xml", "w")
	print >> file, template % s

for s in strain:
        file = open("./homolog/" + s + "_" + s + ".xml", "w")
        print >> file, paralog_template % (s, s)

for i, s1 in enumerate(strain):
	for j, s2 in enumerate(strain[i + 1:]):
		file = open("./homolog/" + s1 + "_" + s2 + ".xml", "w")
		print >> file, hmlg_template % (s1, s2)
