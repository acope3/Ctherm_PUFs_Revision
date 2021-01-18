from bs4 import BeautifulSoup
import glob
import re
import os

files = glob.glob("*/report.html")
pat = re.compile("(ADU[0-9]+)|(WP_[0-9]+)")
for f in files:
	print f
	folder = f.split("/")[0]
	protein = pat.search(f).group()
	protein = protein[:-1]+"."+protein[-1]
	counts = {}
	names =set([])
	with open(f) as fin, open(folder+"/"+protein+"_structures.txt",'w') as out:
		soup = BeautifulSoup(fin)
		table = soup.find_all("tr")
		for row in table[1:]:
			cols = row.find_all("td")
			structure = str(cols[9]).lower()[4:-5]
			names.add(structure)
			current = counts.get(structure)
			if current == None:
				counts[structure] = 1
			else:
				counts[structure] = current+1
		out.write("Structure\tNumber\n")
		for i in names:
			out.write("\t".join([i,str(counts.get(i))]))
			out.write("\n")

