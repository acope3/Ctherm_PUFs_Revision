import pandas as pd
from Bio import SeqIO

def parsePANNZER2(target_file):
	df = pd.read_csv(target_file,sep="\t",header=0,dtype={"goid":"string"})
	df["goid"] = df.goid.astype(str)
	print len(set(df["qpid"].values))
	df_grouped = df.loc[:,["qpid","goid"]].groupby(["qpid"])["goid"].apply(lambda x: ",".join("GO:"+x)).reset_index()
	return df_grouped

def parseINTERPRO(target_file):
	df = pd.read_csv(target_file,sep="\t",header=0,dtype={"GO_annot":"string"})
	#df.fillna(value="NA",inplace=True)
	print len(set(df["Protein_Accession"].values))
	df.dropna(axis=0,how="any",subset=["GO_annot"],inplace=True)
	df_grouped = df.loc[:,["Protein_Accession","GO_annot"]].groupby(["Protein_Accession"])["GO_annot"].apply(lambda x: ",".join(x)).reset_index()
	df_grouped["GO_annot"] = df_grouped["GO_annot"].apply(lambda x: ",".join(x.split("|")))
	return df_grouped


def getProtWNoGo():
	prot = SeqIO.parse(open("Data/FASTA/Ctherm_RefSeq_Genbank_updated.fasta"),"fasta")
	go = pd.read_csv("Data/GO/ctherm_protein_to_go_map.tsv",sep="\t",header=0,index_col=0)
	with open("Data/FASTA/ctherm_no_go.fasta",'w') as out:
		for fasta in prot:
			prot_id = fasta.id
			try:
				value = go.loc[prot_id,"GO"]
				if pd.isnull(value):
					out.write(">"+fasta.description+"\n"+str(fasta.seq)+"\n")
			except KeyError:
				out.write(">"+fasta.description+"\n"+str(fasta.seq)+"\n")

#getProtWNoGo()


# pan = parsePANNZER2("Data/PANNZER2/GO.out")
# pan["ProteinID"] = pan.loc[:,"qpid"]
# interpro = parseINTERPRO("Data/GO/Ctherm_Inter_results2.tsv")
# interpro["ProteinID"] = interpro.loc[:,"Protein_Accession"]

# eggnog = pd.read_csv("Data/GO/Ctherm_RefSeq_Genbank_NR_emapper_annotations.tab",sep="\t",header=0,usecols=["ProteinID","GO_terms"])

# tmp = pan.merge(interpro,how="outer",left_on="ProteinID",right_on="ProteinID")
# tmp2 = tmp.merge(eggnog,how="outer",left_on="ProteinID",right_on="ProteinID")
# tmp2 = tmp2.loc[:,["ProteinID","goid","GO_annot","GO_terms"]]

# tmp2.fillna(value="NA",inplace=True)
# tmp2["GO"] = tmp2['goid']+','+tmp2['GO_annot']+','+tmp2['GO_terms']

# tmp2["GO"] = tmp2["GO"].apply(lambda row: ",".join(list(set(row.split(",")))))
# tmp2["GO"] = tmp2["GO"].str.replace(r",*NA,*","")

# tmp2[["ProteinID","GO"]].to_csv("Data/GO/ctherm_protein_to_go_map.tsv",sep="\t",header=True,index=False,quoting=False)

current = pd.read_csv("Data/GO/ctherm_protein_to_go_map.tsv",sep="\t",header=0)
pan = parsePANNZER2("Data/PANNZER2/No_go_reanalysis/GO.out")
pan["ProteinID"] = pan.loc[:,"qpid"]

eggnog = pd.read_csv("Data/GO/No_go/ctherm_no_go.fa.emapper.annotations",sep="\t",header=0)
eggnog["ProteinID"] = eggnog.loc[:,"query_name"]



tmp = pan.merge(eggnog,on="ProteinID",how="outer")
current = current.merge(tmp,on="ProteinID",how="outer")

current.fillna(value="NA",inplace=True)
current["GO"] = current["GO"]+","+current['goid']+','+current['GOs']

current["GO"] = current["GO"].apply(lambda row: ",".join(list(set(row.split(",")))))
current["GO"] = current["GO"].str.replace(r",*NA,*","")


current[["ProteinID","GO"]].to_csv("Data/GO/ctherm_protein_to_go_map_updated.tsv",sep="\t",header=True,index=False,quoting=False)

