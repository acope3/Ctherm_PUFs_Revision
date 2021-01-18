import pandas as pd

def getProtFunctions():
	pufs = []
	function = {}
	with open("../Data/FASTA/Ctherm_RefSeq_Genbank_updated.fasta") as fin:
		for line in fin:
			if line[0] == ">":
				line_spt = line.split()
				if ("hypothetical" in line) or ("DUF" in line):
					pufs.append(line_spt[0][1:])
				ind = line.find("[")
				if ind != -1:
					ind_space = line.find(" ")
					function[line_spt[0][1:]] = line[ind_space+1:ind-1]
				else:
					function[line_spt[0][1:]] = "None"
	return function,pufs


pufs = []
with open("../Data/measured_pufs.txt") as fin:
	for line in fin:
		pufs.append(line.strip())

ortho = pd.read_csv("Results_Jul22_3/Orthogroups/Orthogroups.tsv",sep="\t",header=0,index_col=0)
ctherm = ortho.loc[:,"cthermocellum_protein"]
ctherm.dropna(how="any",inplace=True)
ortho_pufs = set([])
total = 0
for i in pufs:
	group = ctherm[ctherm.str.contains(i)]
	if len(group)>0:
		if "," in group.values[0]:
			print group.index.values[0]
			total+=1
			ortho_pufs.add(group.index.values[0])
print total,len(ortho_pufs),len(pufs)
function,pufs = getProtFunctions()
for i in ortho_pufs:
	x = ctherm.loc[i]
	print i, x
	prot = x.split(",")
	for j in prot:
		print function.get(j.strip())