template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "%s_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "strand" />
	</Dataset>
</Query>'''

for line in open("mice.txt"):
	strain = line.strip().split()[1].split('_')[0]
	file = open("./gene/" + strain + "_gene.xml", "w")
	print >> file, template % strain
