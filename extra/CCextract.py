import pandas
import numpy as np

my_list = np.zeros(10000000, dtype=bool)
with open("CCnode.txt", "r") as f:
	for line in f:
		my_list[int(line.strip().split("\t")[0])] = True

with open("refined_nodes.tsv", "w") as wf:
	df = pandas.read_csv("/networks/net_refined_level1_thr3_20190214.txt", sep='\t', header=None)
	for i, row in enumerate(df.values):
		a, b = row[0], row[1]
		a, b = int(a), int(b)
		if my_list[a] and my_list[b]:
			wf.write("%d\t%d\n"%(a, b))
			print(a, b, "  >>> added")
