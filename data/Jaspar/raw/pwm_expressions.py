import csv
import numpy as np

# JASPAR pwms
with open(r'.\JASPAR_CORE_nonredundant_pfm_vertebrates.txt', 'r') as f:
    pwms_raw = f.read().lower()

# TFs ids to their names
with open(r'.\Ensembl_v65.Gencode_v10.ENSG.gene_info') as f:
    reader = csv.reader(f, delimiter='\t')
    names_rows = [row for row in reader]
    names_dict = dict([[row[0], row[6]] for row in names_rows])

# expression of TFs, by id
with open(r'.\57epigenomes.N.pc') as f:
    reader = csv.reader(f, delimiter='\t')
    expressions_rows = [row for row in reader]

header = expressions_rows[0]
data_rows = np.array(expressions_rows[1:])

names = np.vectorize(lambda name: names_dict[name])(data_rows[:, :1])
names_mask = [name.lower() in pwms_raw and name not in ['', 'NA'] for name in names[:, 0]]
names = names[names_mask, :1]
data = data_rows[names_mask, 1:-1].astype(float)
data_norm = data / data.sum(axis=1)[:, np.newaxis]
data_str = np.round(data_norm, 3).astype(str)
output_data = np.vstack([header, np.hstack([names, data_str])])


# write TFs expression by name
np.savetxt(r'.\TF_expressions.txt', output_data, delimiter='\t', fmt='%s')
