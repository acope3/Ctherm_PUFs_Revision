from Bio import SeqIO
import re

def updateOldAnnotations():
	old_proteome =  SeqIO.to_dict(SeqIO.parse(open("FASTA/Ctherm_RefSeq_Genbank_NR_cRAP.fasta"),"fasta"))

	new_refseq =  SeqIO.to_dict(SeqIO.parse(open("FASTA/ctherm_refseq_protein.faa"),"fasta"))
	new_gen = SeqIO.to_dict(SeqIO.parse(open("FASTA/ctherm_GenBank_protein.faa"),"fasta"))

	old_ids = old_proteome.keys()


	with open("Ctherm_RefSeq_Genbank_updated.fasta",'w') as out:
		for prot in old_ids:
			try: 
				if prot[0:2] == "AD":
					prot_info = new_gen.get(prot)
					out.write(">"+prot_info.description+"\n")
					out.write(str(prot_info.seq)+"\n")
				elif prot[0:2] == "WP":
					prot_info = new_refseq.get(prot)
					out.write(">"+prot_info.description+"\n")
					out.write(str(prot_info.seq)+"\n")
				else:
					print prot
			except AttributeError:
				prot_info = old_proteome.get(prot)
				out.write(">"+prot_info.description+"\n")
				out.write(str(prot_info.seq)+"\n")


def getProtID(rec,pat):
	value = pat.search(rec.description)
	if value:
		return value.group(1)
	else:
		return rec.id

def filterPseudo(fna,target):
	cds = SeqIO.parse(open(fna),"fasta")
	with open(target,'w') as out:
		for fasta in cds:
			if ("[pseudo=true]" not in fasta.description) and ("transposase" not in fasta.description):
				out.write(">"+fasta.description+"\n")
				out.write(str(fasta.seq)+"\n")


def createCDSFile():
	old_proteome =  SeqIO.to_dict(SeqIO.parse(open("FASTA/Ctherm_RefSeq_Genbank_NR_cRAP.fasta"),"fasta"))

	pat = re.compile(r"\[protein_id=([A-Z\.0-9_]+)\]")

	new_refseq =  SeqIO.parse(open("FASTA/ctherm_refseq_cds_from_genomic_remove_pseudo_transposase.fna"),"fasta")
	new_gen = SeqIO.parse(open("FASTA/ctherm_GenBank_cds_from_genomic_remove_pseudo_transposase.fna"),"fasta")

	with open("Ctherm_RefSeq_Genbank_codons.fasta",'w') as out:
		for fasta in new_refseq:
			value = pat.search(fasta.description)
			if value:
				if value.group(1) in old_proteome.keys():
					out.write(">"+fasta.description+"\n")
					out.write(str(fasta.seq)+"\n")
		for fasta in new_gen:
			value = pat.search(fasta.description)
			if value:
				if value.group(1) in old_proteome.keys():
					out.write(">"+fasta.description+"\n")
					out.write(str(fasta.seq)+"\n")

def pullRiboSeq():
	cds = SeqIO.parse(open("FASTA/Ctherm_RefSeq_Genbank_codons.fasta"),"fasta")
	with open("FASTA/ribo_prot.fasta",'w') as out:
		for fasta in cds:
			if "ribosomal" in fasta.description:
				out.write(">"+fasta.description+"\n")
				out.write(str(fasta.seq)+"\n")




#filterPseudo("FASTA/ctherm_refseq_cds_from_genomic.fna","FASTA/ctherm_refseq_cds_from_genomic_remove_pseudo_transposase.fna")
#filterPseudo("FASTA/ctherm_GenBank_cds_from_genomic.fna","FASTA/ctherm_GenBank_cds_from_genomic_remove_pseudo_transposase.fna")
#createCDSFile()
pullRiboSeq()