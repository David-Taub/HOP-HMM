import csv
with open(r'.\Ensembl_v65.Gencode_v10.ENSG.gene_info') as f:
    reader = csv.reader(f, delimiter='\t')
    names_rows = [row for row in reader]
    names_dict = dict([[row[0], row[6]] for row in names_rows])


with open(r'.\57epigenomes.N.pc') as f:
    reader = csv.reader(f, delimiter='\t')
    expressions_rows = [row for row in reader]
    expressions_rows_named = [[names_dict[row[0]]] + row[1:] for row in expressions_rows[1:]]
    expressions_rows_named.insert(0, expressions_rows[0])
    expressions_rows_named.insert(0, expressions_rows[0])

with open(r'.\TF_expressions.txt', 'w') as f:
    [f.write('\t'.join(row) + '\n') for row in expressions_rows_named]
