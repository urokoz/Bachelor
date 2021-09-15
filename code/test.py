infile = open("Data/25885_NetMHCIIpan.xls","r")

infile.readline()
infile.readline()
pep_HLA_dict = dict()
old_pep = 0
for i, line in enumerate(infile):
    line = line.split()

    cur_pep = line[2]
    if old_pep != cur_pep:
        if old_pep:
            pep_HLA_dict[old_pep] = [2 if a<1 else 1 if a<5 else 0 for a in pep_HLA_dict[old_pep]]

        old_pep = cur_pep

    HLA_bind_rank = [float(line[i]) for i in range(6,25,3)]

    if cur_pep in pep_HLA_dict:
        pep_HLA_dict[cur_pep] = [min(old, new) for old, new in zip(pep_HLA_dict[cur_pep], HLA_bind_rank)]
    else:
        pep_HLA_dict[cur_pep] = HLA_bind_rank

pep_HLA_dict[old_pep] = [2 if a<1 else 1 if a<5 else 0 for a in pep_HLA_dict[old_pep]]
