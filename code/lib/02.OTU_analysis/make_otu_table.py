import sys


###<otu><tax><outfile>
(otup,tax,outfile)=sys.argv[1:]
res = {}
with open(tax,'r') as taxfile:
	for line in taxfile:
		if not line or line.startswith('#'):
			continue
		line = line.rstrip("\n")
		fields = line.split('\t')
		otu = fields[0].split(' ')[0]
		res[otu] = fields[1].replace(';','; ')
out=open(outfile,'w')
with open(otup,'r') as otufile:
	line=otufile.readline()
	line='%s\ttaxonomy\n' %line.strip()
	out.write(line)
	for line in otufile:
		if not line or line.startswith('#'):
                        continue
		line=line.rstrip("\n")
		fields = line.split('\t')
		line='%s\t%s\n' %(line,res[fields[0]])
		out.write(line)
out.close()
