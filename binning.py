# module load python/3.8.1
# run with `python3`, not `python`
# also needs seqkit loaded to pull out the subsetted sequences
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns

markdir = "/home/CAM/mstukel/AHE/Kikihia/blastdb"

#single file test version
testfile = "I6899contigs.fasta"
testpath = os.path.join(markdir, testfile)

pd_data = {"assembly": [], "node": [], "length": [], "coverage": []}

for record in SeqIO.parse(testpath, "fasta"):
    # seq name format is NODE_2876_length_652_cov_2.152429
    x = record.id.split("_")
    pd_data["assembly"].append(testfile)
    pd_data["node"].append(x[1])
    pd_data["length"].append(int(x[3]))
    pd_data["coverage"].append(float(x[-1]))

df = pd.DataFrame(data = pd_data)

# df.astype({"length": "int32", "coverage": "float"})

print(df)

#example matplotlib commands
#in the future maybe try binning with GC content?

#lengths histogram
df["length"].plot(kind="hist")
plt.yscale('log', nonpositive='clip')  #log scale y
plt.savefig("./testfig_length.png")
plt.clf()   # clear figure

#coverage histogram
df["coverage"].plot(kind="hist")
plt.yscale('log', nonpositive='clip')  #log scale y
plt.savefig("./testfig_coverage.png")
plt.clf()   # clear figure

#len x coverage scatter
df.plot(x="length",y="coverage", kind="scatter")
plt.savefig("./testfig_lenxcoverage.png")
plt.clf()   # clear figure

#len x coverage scatter, log scaled, transparent points
df.plot.scatter(x="length", y="coverage", alpha = 0.1)
plt.yscale('log', nonpositive='clip')  #log scale y
plt.xscale('log', nonpositive='clip')  #log scale x
plt.savefig("./testfig_lenxcoverage_logscaled.png")
plt.clf()   # clear figure

#### 12/17/20 Five file test version
markdir = "/home/CAM/mstukel/AHE/Kikihia/blastdb"
blastresdir = "/home/CAM/mstukel/AHE/fullAHEloci/rmtaxaout/regrouped/blastresults"

os.listdir(markdir)
assemblyfiles = ["I6885contigs.fasta", "I27905contigs.fasta", "I6899contigs.fasta", "I6843contigs.fasta", "I6916contigs.fasta"]
assemblypaths = [os.path.join(markdir, x) for x in assemblyfiles]

pd_data = {"assembly": [], "node": [], "length": [], "coverage": []}
for path in assemblypaths:
    for record in SeqIO.parse(path, "fasta"):
        # seq name format is NODE_2876_length_652_cov_2.152429
        x = record.id.split("_")
        pd_data["assembly"].append(path.rsplit("/", 1)[1])
        pd_data["node"].append(int(x[1]))
        pd_data["length"].append(int(x[3]))
        pd_data["coverage"].append(float(x[-1]))

df = pd.DataFrame(data = pd_data)
# df = df.set_index("assembly")
del pd_data

blastrespaths = [os.path.join(blastresdir, x) + ".blast" for x in assemblyfiles]

pd2_data = {"assembly": [], "locus": [], "node": [], "length": [], "coverage": []}
for file in blastrespaths:
    result = subprocess.run(args = "sort -u -k1,1 {}".format(file), shell = True, capture_output=True)
    blastfile = result.stdout.decode("utf-8")
    for line in blastfile.splitlines():
        x = line.split("\t")
        nodedata = x[1].split("_")
        pd2_data["assembly"].append(file.rsplit("/", 1)[1][:-6])
        pd2_data["node"].append(int(nodedata[1]))
        pd2_data["length"].append(int(nodedata[3]))
        pd2_data["coverage"].append(float(nodedata[-1]))
        pd2_data["locus"].append(x[0].lstrip("./")[:-4])

best_hits = pd.DataFrame(data = pd2_data)
del pd2_data

df2 = pd.merge(df, best_hits, on=["assembly", "node"], how = "left", indicator = True, suffixes=("", "_y"))
df2["target"] = np.where(df2["_merge"] == "both", True, False)
df2.drop(list(df2.filter(regex = "_")), axis = 1, inplace = True)
df2 = df2.set_index("assembly")

#GRID WITH THE ALIBA-PULLED SEQS HIGHLIGHTED FOR 5 ASSEMBLIES
df2.reset_index(level=0, inplace = True)
grid = sns.FacetGrid(df2, col = "assembly", col_wrap = 3)
grid.map_dataframe(sns.scatterplot, x = "length", y = "coverage", hue = "target", palette = "tab10", alpha = 0.9, s = .5)
grid.set(yscale = "log", xscale = "log")
grid.set_axis_labels("log(length)", "log(coverage)")
grid.set_titles(col_template="{col_name}")
plt.tight_layout()
plt.savefig("./testfig_grid.png", dpi = 300)
plt.clf()

try:  ## you only need to run this if you ran the block above, otherwise it'll throw an error because the index is already set
    df2 = df2.set_index("assembly")
except Exception as e:
    pass

## subsetting parameters here
subset = df2[(df2["length"] > 200) & (df2["coverage"] > 9)]
np.where(df2["target"] == True)[0].size

#reduced dataset size by 97%
# 1-subset.loc["I6885contigs.fasta"].shape[0]/df2.loc["I6885contigs.fasta"].shape[0]
1-subset.shape[0]/df2.shape[0]

#retained 97% of target sequences
np.where(subset["target"] == True)[0].size/np.where(df2["target"] == True)[0].size

#20% of remaining sequences are target sequences
np.where(subset["target"] == True)[0].size/subset.shape[0]
np.where(df2["target"] == True)[0].size/df2.shape[0]


df2.index.unique().tolist()

# makes new files with the subsetted sequences from above
os.makedirs("subsets")
for x in subset.index.unique().tolist():
    working_subset = subset.loc[x]
    # working_subset = working_subset[working_subset["target"] == True]
    names = working_subset.apply(lambda row: "NODE_{n}_length_{l}_cov_{c}".format(n = row["node"], l = row["length"], c = row["coverage"]), axis = 1).tolist()
    newfile = "./subsets/{}_subsetnames.txt".format(x)
    with open(newfile, "w") as f:
        f.write("\n".join(names))
    subprocess.run(args = "seqkit grep -f {n} {m} > {f}".format(n = newfile, m = os.path.join(markdir,x), f = "./subsets/{}_subset.fasta".format(x[:-6])), shell = True, capture_output=False)
