import pandas as pd
import numpy as np
import scipy as sp
import glob
import re
import os.path

def pullDoorDBAnnot(filename,annot_map):
	door_annot = {}
	with open(filename) as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split("\t")
			try:
				ids = annot_map.get(line_spt[2])
				for i in ids:
					door_annot[i] = line_spt[-1]
			except TypeError:
				continue
	return door_annot


def oldToNewAnnot(refseq,genbank):
	annot_map = {}
	with open(refseq) as fin:
		fin.readline()
		pat = re.compile("old_locus_tag=")
		for line in fin:
			line_spt = line.strip().split("\t")
			if line_spt[0] == "gene":
				if pat.search(line_spt[-1]) != None:
					locus = line_spt[-1].split("=")[1]
				else:
					locus = None
			elif line_spt[0] == "CDS" and locus != None:
				annot_map[locus] = [line_spt[10]]
			else:
				continue
	with open(genbank) as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split("\t")
			if line_spt[0] == "gene":
				locus = line_spt[16]
			elif line_spt[0] == "CDS":
				try:
					annot_map[locus].append(line_spt[10])
				except KeyError:
					annot_map[locus] = [line_spt[10]]
			else:
				continue
	return annot_map


def readInterPro(filename,desc_index=12,go_index=13):
	ip_desc = {}
	ip_go = {}
	with open(filename) as fin:
		fin.readline()
		for line in fin:
			line_spt = line.split("\t")
			prot = line_spt[0]
			prot_desc = ip_desc.get(prot)
			prot_go = ip_go.get(prot)
			if line_spt[desc_index] != '':
				if prot_desc != None:
					prot_desc.add(line_spt[desc_index].strip())
				else:
					ip_desc[prot] = set([line_spt[desc_index].strip()])
			go_terms = line_spt[go_index].split("|")
			if go_terms != ['']:
				if prot_go != None:
					for go in go_terms:
						prot_go.add(go.strip())
				else:
					ip_go[prot] = set(go_terms)
	return ip_desc,ip_go



def getProtFunctions():
	pufs = []
	function = {}
	with open("Data/FASTA/Ctherm_RefSeq_Genbank_updated.fasta") as fin:
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


def parsePhylogeneticTrees():
	files = glob.glob("Evolutionary/Phylogenetic_trees/Ctherm_PUFs_allTrees/*")
	tree_functions = {}
	tip_label_pat = re.compile(r"tr\|[A-Z0-9]+\|[A-Z0-9]+_[A-Z0-9]+_([A-Za-z0-9_]+):[0-9]\.[0-9]+")
	for f in files:
		with open(f) as fin:
			prot = f.split("_MAFFT_crude.fasta_trim_tree")[0].split("/")[-1]
			tree = fin.readline()
			tip_label = tip_label_pat.finditer(tree)
			functions = set([])
			for i in tip_label:
				functions.add(i.group(1))
			tree_functions[prot] = functions
	return tree_functions

def getOldLocusTag(refseq,genbank):
	annot_map = {}
	with open(refseq) as fin:
		fin.readline()
		pat = re.compile("old_locus_tag=")
		for line in fin:
			line_spt = line.strip().split("\t")
			if line_spt[0] == "gene":
				if pat.search(line_spt[-1]) != None:
					locus = line_spt[-1].split("=")[1]
				else:
					locus = None
			elif line_spt[0] == "CDS" and locus != None:
				annot_map[line_spt[10]] = locus
			else:
				continue
	with open(genbank) as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split("\t")
			if line_spt[0] == "gene":
				locus = line_spt[16]
			elif line_spt[0] == "CDS":
				annot_map[line_spt[10]] = locus
			else:
				continue
	return annot_map



def parseSwissModel():
	files = glob.glob("SwissModel_structures/*/*_structures.txt")
	structures = {}
	for f in files:
		protein = f.split("/")[-1].split("_structures.txt")[0]
		with open(f) as fin:
			desc = []
			fin.readline()
			for line in fin:
				line_spt = line.strip().split("\t")
				desc.append((line_spt[0],line_spt[1]))
			structures[protein] = desc
	return structures


def pullHMMAnnotation(filename):
	df = pd.read_table(filename,sep="\t",header=0,index_col=0)
	return df.iloc[:,-1].T.to_dict()

def readInterPro(filename,desc_index=12):
	ip_desc = {}
	with open(filename) as fin:
		fin.readline()
		for line in fin:
			line_spt = line.split("\t")
			prot = line_spt[0]
			prot_desc = ip_desc.get(prot)
			if line_spt[desc_index].strip() != '':
				if prot_desc != None:
					prot_desc.add(line_spt[desc_index].strip())
				else:
					ip_desc[prot] = set([line_spt[desc_index].strip()])
		
	return ip_desc



def readSignalPOutput(filename):
	signalp = []
	with open(filename) as fin:
		for line in fin:
			signalp.append(line.split()[0])
	return signalp

def readTMHMMOutput(filename):
	tm = {}
	with open(filename) as fin:
		for line in fin:
			line_spt = line.strip().split("\t")
			tm[line_spt[0]] = line_spt[4].split("=")[1]
	return tm




def findPUFinClust(puf,clust_file):
	with open(clust_file) as fin:
		for line in fin:
			line_spt = line.split()
			try: 
				cluster_num = line_spt.index(puf)
				return cluster_num
			except ValueError:
				continue
	return -1

def getGOEnrichment(go_file):
	df = pd.read_csv(go_file,header=0,sep="\t")
	df_sig = df.loc[df.weight01Fisher < 0.05,"GO.ID"].values
	return df_sig

def getKEGGEnrichment(kegg_file):
	df = pd.read_csv(kegg_file,header=0,sep="\t")
	df_sig = df.loc[df.qvalue < 0.05,"Description"].values
	return df_sig

def getDifferentiallyExpressed(diff_exp_file):
	df = pd.read_csv(diff_exp_file,header=0,sep="\t",index_col=0)
	df_sig = df.loc[(df.loc[:,"adj.P.Val"] < 0.05) & (df.PUF == "Yes"),:]
	return df_sig


def createPUFsTable():
	function,pufs = getProtFunctions()

	genome_sequences = function.keys()
	old_locus_ids = getOldLocusTag("Data/GCF_000184925.1_ASM18492v1_feature_table.txt","Data/GCA_000184925.1_ASM18492v1_feature_table.txt")
	new_annot = oldToNewAnnot("Data/GCF_000184925.1_ASM18492v1_feature_table.txt","Data/GCA_000184925.1_ASM18492v1_feature_table.txt")
	door_annot = pullDoorDBAnnot("Data/DOOR/Ctherm_1313_operon.tab",new_annot)
	structures = parseSwissModel()
	trees = parsePhylogeneticTrees()
	hmm = pullHMMAnnotation("Data/GO/Ctherm_RefSeq_Genbank_NR_emapper_annotations.tab")
	interpro_desc = readInterPro("Data/GO/Ctherm_Inter_results2.tsv",desc_index=12)
	interpro_pathway = readInterPro("Data/GO/Ctherm_Inter_results2.tsv",desc_index=14)
	pannzer_annot = pd.read_csv("Data/puf_duf_new_annotations.txt",sep="\t",header=0,index_col=0)
	kegg = pd.read_csv("Data/GO/Ctherm_PUF_BlastKOALA.txt",sep="\t",header=0,index_col=0)
	ec = pd.read_csv("Data/ec_pufs.txt",sep="\t",header=0,index_col=0)
	ec_function = pd.read_csv("Data/ec_to_function.txt",sep="\t",header=0,index_col=0)
	protein_go = pd.read_csv("Data/GO/ctherm_protein_to_go_map_updated.tsv",sep="\t",header=0,index_col=0)
	go_terms = pd.read_csv("Data/GO/go_terms.tsv",sep="\t",header=0,index_col=0)
	sp_prot = readSignalPOutput("Data/sp.gff")
	tm_prot = readTMHMMOutput("Data/tmhmm_results.txt")
	operons = pd.read_csv("Data/DOOR/operon_member_table.txt",header=0,sep="\t")
	
	diff_exp_t1 = getDifferentiallyExpressed("DiffExp/Tables/hpt_vs_LL1210_early_log.tsv")
	diff_exp_t2 = getDifferentiallyExpressed("DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv")
	diff_exp_t3 = getDifferentiallyExpressed("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv")

	diff_coexpress = pd.read_csv("DiffCoExpress/differential_coexpression_strains.tsv",sep="\t",header=0)
	diff_coexpress_sig = diff_coexpress.loc[(diff_coexpress.loc[:,"q.value"] < 0.05) & (diff_coexpress.PUF == "Yes"),"Protein"].values

	coevolution = pd.read_csv("Evolutionary/2020-07-22-phylo_prof_significant_pufs_with_prot_ids.tsv",sep="\t",header=0)

	

	df = pd.DataFrame(columns=["Old Locus ID",
		"Cluster",
		"Operon Proteins",
		"Operon Proteins (Old IDs)" ,
		"Operon Functions",
		"Signal Peptide",
		"Transmembrane Region",
		"eggNOG Annotation",
		"PANNZER2 Annotation",
		"InterProScan Annotation",
		"InterProScan Pathway",
		"BlastKOALA",
		"EC Number",
		"EC Function",
		"GO_terms",
		"LogFold_Change_Early_Log",
		"LogFold_Change_Mid_Log",
		"LogFold_Change_Late_Log",
		"Differentially_coexpressed",
		"hpt_Clust_GO_Biological_Process",
		"hpt_Clust_GO_Molecular Function",
		"LL1210_Clust_GO_Biological_Process",
		"LL1210_Clust_GO_Molecular Function",
		"Structures",
		"Coevolving",
		"Gene_trees"])
	measured_pufs = pd.read_csv("Data/measured_pufs.txt",delim_whitespace=True,header=None)
	for protein in measured_pufs[[0]].loc[:,0].values:
		


		hpt_bp_enrich = "NA"
		hpt_mf_enrich = "NA"
		hpt_kegg_enrich = "NA"
		hpt_kegg_mod_enrich = "NA"
		ll_bp_enrich = "NA"
		ll_mf_enrich = "NA"
		ll_kegg_enrich = "NA"
		ll_kegg_mod_enrich = "NA"


		clust_index_hpt = str(findPUFinClust(protein,"Clust/hpt/Clusters_Objects.tsv"))
		clust_index_ll = str(findPUFinClust(protein,"Clust/LL1210/Clusters_Objects.tsv"))
		
		if clust_index_hpt != "-1":
			
			try:
				clust_go = getGOEnrichment("Clust/hpt/GO/go_enrichment_BP_cluster_"+clust_index_hpt+".tsv")
				if len(clust_go) > 0:
					hpt_bp_enrich = ";".join(clust_go)
			except IOError:
				print "File not found"
			try:	
				clust_go = getGOEnrichment("Clust/hpt/GO/go_enrichment_MF_cluster_"+clust_index_hpt+".tsv")
	
				if len(clust_go) > 0:
					hpt_mf_enrich = ";".join(clust_go)
			except IOError:
				print "File not found"	
			try:
				clust_kegg = getKEGGEnrichment("Clust/hpt/KEGG/kegg_enrichment_clust_"+clust_index_hpt+".tsv")

				if len(clust_kegg) > 0:
					hpt_kegg_enrich = ";".join(clust_kegg)
			except IOError:
				print "File not found"
			try:
				clust_kegg = getKEGGEnrichment("Clust/hpt/KEGG/kegg_enrichment_module_clust_"+clust_index_hpt+".tsv")
				if len(clust_kegg) > 0:
					hpt_kegg_mod_enrich = ";".join(clust_kegg)
			
			except IOError:
				print "File not found"
			
		else:
			clust_index_hpt = "NA"
		

		if clust_index_ll != "-1":
			try:
				clust_go = getGOEnrichment("Clust/LL1210/GO/go_enrichment_BP_cluster_"+clust_index_ll+".tsv")
				if len(clust_go) > 0:
					ll_bp_enrich = ";".join(clust_go)
			except IOError:
				print "File not found"
			try:
				clust_go = getGOEnrichment("Clust/LL1210/GO/go_enrichment_MF_cluster_"+clust_index_ll+".tsv")
				if len(clust_go) > 0:
					ll_mf_enrich = ";".join(clust_go)
			except IOError:
				print "File not found"
			try:
				clust_kegg = getKEGGEnrichment("Clust/LL1210/KEGG/kegg_enrichment_clust_"+clust_index_ll+".tsv")
				if len(clust_kegg) > 0:
					ll_kegg_enrich = ";".join(clust_kegg)
			except IOError:
				print "File not found"
			try:
				clust_kegg = getKEGGEnrichment("Clust/LL1210/KEGG/kegg_enrichment_module_clust_"+clust_index_ll+".tsv")	
				if len(clust_kegg) > 0:
					ll_kegg_mod_enrich = ";".join(clust_kegg)
			except IOError:
				print "File not found"
		
		else:
			clust_index_ll = "NA"

		

		clust = str(clust_index_hpt) + ";" + str(clust_index_ll)
		protein_old_id = old_locus_ids.get(protein)
		puf_operon = []
		operon_members = []
		old_ids = []
		for k,row in operons.iterrows():
			if protein in operons.iloc[k,1].split(";"):
				for i in operons.iloc[k,1].split(";"):
					if i in genome_sequences:
						puf_operon.append(door_annot.get(i))
						operon_members.append(i)
						old_ids.append(old_locus_ids.get(i))
				break
		puf_operon =";".join(puf_operon)
		operon_members = ";".join(operon_members)
		old_ids = ";".join(old_ids)
		hmm_annot = str(hmm.get(protein))
		try:
			puf_pannzer = pannzer_annot.loc[protein,"Pannzer"]
		except KeyError:
			puf_pannzer = "NA"
		
		puf_interpro_desc = interpro_desc.get(protein)
		if puf_interpro_desc != None:
			puf_interpro_desc = ";".join(list(puf_interpro_desc))
		else:
			puf_interpro_desc = "NA"
		puf_interpro_pathway = interpro_pathway.get(protein)
		if puf_interpro_pathway != None:
			puf_interpro_pathway = ";".join(list(puf_interpro_pathway))
		else:
			puf_interpro_pathway = "NA"
		
		try:
			puf_ec = ec.loc[protein,"EC"]

		except KeyError:
			puf_ec = "NA"
		
		try:
			puf_kegg = kegg.loc[protein,"KEGG_Term"]
		except KeyError:
			puf_kegg = "NA"

		if protein in sp_prot:
			sp = True
		else:
			sp = False
		tm = tm_prot.get(protein)
	
		try:
			el_fc = diff_exp_t1.loc[protein,"logFC"]
		except KeyError:
			el_fc = "NA"
		try:
			ml_fc = diff_exp_t2.loc[protein,"logFC"]
		except KeyError:
			ml_fc = "NA"
		try:
			ll_fc = diff_exp_t3.loc[protein,"logFC"]
		except KeyError:
			ll_fc = "NA"
		
		if protein in diff_coexpress_sig:
			diff_co = "True"
		else:
			diff_co = "No"

		
		
		coevol_puf = coevolution.loc[coevolution.PUF_Protein_Ids == protein,"Target_Protein_Ids"].values
		coevol_puf = map(lambda x: x.split(", "),coevol_puf)
		coevol_puf = [item for sublist in coevol_puf for item in sublist]
		
		coevol_puf = map(lambda x: str(function.get(x)),coevol_puf)
		coevol_puf = ";".join(list(coevol_puf))

		puf_tree = trees.get(protein)
		if puf_tree != None:
			puf_tree = ";".join(list(puf_tree))
		else:
			puf_tree = "NA"
		try:
			go_puf = str(protein_go.loc[protein,"GO"])
			if go_puf != "nan":
				go_puf_desc = []
				for i in go_puf.split(","):
					go_puf_desc.append(go_terms.loc[i,"Name"])
				go_puf_desc = ";".join(go_puf_desc)
			else:
				go_puf_desc = "NA"
		except KeyError:
			go_puf_desc = "NA"
		models = structures.get(protein)
		tmp = []
		if models != None:
			for m in models:
				tmp.append(m[0]+"("+m[1]+")")
			models = ";".join(tmp)
		else:
			models = "NA"
		value  = pd.DataFrame({
			"Cluster":clust,
			"Operon Proteins":operon_members,
			"Operon Functions":puf_operon,
			"Signal Peptide":sp,
			"Transmembrane Region":tm,
			"eggNOG Annotation":hmm_annot,
			"PANNZER2 Annotation":puf_pannzer,
			"InterProScan Annotation":puf_interpro_desc,
			"InterProScan Pathway":puf_interpro_pathway,
			"BlastKOALA":puf_kegg,
			"EC Number":puf_ec,
			"GO_terms":go_puf_desc,
			"LogFold_Change_Early_Log":el_fc,
			"LogFold_Change_Mid_Log":ml_fc,
			"LogFold_Change_Late_Log":ll_fc,
			"Differentially_coexpressed":diff_co,
			"hpt_Clust_GO_Biological_Process":hpt_bp_enrich,
			"hpt_Clust_GO_Molecular Function":hpt_mf_enrich,
			"LL1210_Clust_GO_Biological_Process":ll_bp_enrich,
			"LL1210_Clust_GO_Molecular Function":ll_mf_enrich,
			"Structures":models,
			"Coevolving":coevol_puf,
			"Gene_trees":puf_tree},index=[protein])
		df = df.append(value)
	df
	return df

pufs_table = createPUFsTable()
pufs_table.to_csv("pufs_table_11_25_2020.tsv",sep="\t",header=True,quoting=False,index=True,
	columns=["eggNOG Annotation",
		"PANNZER2 Annotation",
		"InterProScan Annotation",
		"InterProScan Pathway",
		"BlastKOALA",
		"EC Number",
		"EC Function",
		"GO_terms",
		"Operon Proteins",
		"Operon Functions",
		"LogFold_Change_Early_Log",
		"LogFold_Change_Mid_Log",
		"LogFold_Change_Late_Log",
		"Differentially_coexpressed",
		"Cluster",
		"hpt_Clust_GO_Biological_Process",
		"hpt_Clust_GO_Molecular Function",
		"LL1210_Clust_GO_Biological_Process",
		"LL1210_Clust_GO_Molecular Function",
		"Signal Peptide",
		"Transmembrane Region",
		"Structures",
		"Coevolving",
		"Gene_trees"])
