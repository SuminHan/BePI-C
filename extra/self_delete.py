import pandas

df = pandas.read_csv("/networks/WikiLink.tsv", sep='\t', header=None)

with open("WikiLink_refined.tsv", "w") as wf:
	print("Opened File...")
	for i, row in enumerate(df.values):
		a, b = int(row[0]), int(row[1])
		if a != b:
			wf.write("%d\t%d\n"%(a,b))
		else:
			print("%d deleted"%i)
