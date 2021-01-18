
def parsePannzer():
	with open("No_go_reanalysis/anno.out") as fin,open("../puf_duf_new_annotations.txt",'a') as out,open("../ec_pufs.txt",'a') as out_2:
		fin.readline()
		new_annot = {}
		ec_annot = {}
		for line in fin:
			line_spt = line.strip().split("\t")
			if line_spt[1] == "original_DE":
				current_gene = line_spt[0]
			#elif line_spt[0] == current_gene and line_spt[1] != "qseq" and line_spt[1] != "EC_ARGOT":
			elif line_spt[0] == current_gene and line_spt[1] == "DE":
				annot = new_annot.get(current_gene)
				if annot == None:
					new_annot[current_gene] = [(line_spt[5],line_spt[3])]
				else:
					annot.append((line_spt[5],line_spt[3]))
			elif line_spt[0] == current_gene and line_spt[1] == "EC_ARGOT":
				annot = ec_annot.get(current_gene)
				if annot == None:
					ec_annot[current_gene] = [line_spt[4]]
				else:
					annot.append(line_spt[4])
		print ec_annot
		for key in new_annot.keys():
			out.write(key+"\t")
			tmp = new_annot.get(key)
			funcs = []
			for j in tmp:
				funcs.append(j[0] + "("+j[1]+")")
			out.write(";".join(funcs))
			out.write("\n")
		for key in ec_annot.keys():
			out_2.write(key+"\t")
			tmp = ec_annot.get(key)
			funcs = set([])
			for j in tmp:
				print j
				funcs.add(j.strip())
			out_2.write(";".join(list(funcs)))
			out_2.write("\n")

def parseEC():
	with open("../enzyme.dat") as fin,open("../Data/ec_to_function.txt",'w') as out:
		out.write("EC\tFunction\n")
		for line in fin:
			line_spt = line.split("   ")
			if line_spt[0] == "ID":
				out.write(line_spt[1].strip())
			if line_spt[0] == "DE":
				out.write("\t"+line_spt[1])
parsePannzer()
