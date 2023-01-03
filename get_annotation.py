import pandas as pd
import sys

annotation = pd.read_csv(sys.argv[1], sep = '\t', usecols = [0, 1, 2, 3, 4, 6, 8], skiprows = 5, \
                        names = ['chr', 'source', 'type', 'start', 'end', 'strand', 'attribute'], dtype = {'chr': str})
annotation = annotation.dropna(subset = ['attribute'])

transcript = sys.argv[2]

match = annotation[annotation['attribute'].str.isin(transcript)]
print(match)
