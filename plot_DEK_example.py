import ArcPlot
import RNAStructureObjects 

# initiate object
plot = ArcPlot.ArcPlot()

# add sequence
plot.addFasta('ExampleData/DEKexample/DEK_seq.fa')

# load in pairing probabilities and add them to lower panel
pair_probs = RNAStructureObjects.DotPlot('ExampleData/DEKexample/DEK_invivo_dotplot.dp')
plot.addPairProb(pair_probs, panel=-1)

# generate Matplotlib fig and axis objects
fig, axT, axB = plot.writePlot(write=False)

# read in the phylop data
phylodata = []
with open('ExampleData/DEKexample/dek_3utr_phylop_parsed.txt') as inp:
    for line in inp:
        phylodata.append(float(line.split()[1]))

# plot phylo data in top panel using Matplotlib fill_between function
axT.fill_between(range(len(phylodata)), 0, phylodata)

# format axT using the Matplotlib API 
axT.axes.get_yaxis().set_visible(True)
axT.set_ylabel('PhyloP')
axT.set_ylim(-2,15)

# save the figure
fig.savefig('dek_phylop_invivo.pdf')

