import sys
import pysam
import re

def pull_xml_record(x):
    return(x.split(">")[1].split("<")[0])


seq_ids = {}

with pysam.FastxFile(sys.argv[1]) as fh:
    for entry in fh:
        try:
            if 'locus_tag' in entry.comment:
                comment = entry.comment.split()
                locus_tag =  [x for x in comment if 'locus_tag' in x][0]
                locus_tag = re.sub(pattern = "\[locus_tag=", repl = "", string = locus_tag)
                locus_tag = re.sub("\]", "", locus_tag)
                seq_ids[locus_tag] = entry.name
            else:
                continue
        except:
            continue
        


genus, species, geneid, lineage, geneid_flag, gene_tag = '', '', '', '', 0, ''
with open(sys.argv[2], 'r') as fh:
    for line in fh:
        line = line.strip()
        if '<BioSource>' in line:
            if gene_tag != '':
                print(",".join([seq_ids[gene_tag], gene_tag, genus, species, lineage]))
            genus, species, geneid, lineage, geneid_flag, gene_tag = '', '', '', '', 0, ''
        elif 'BinomialOrgName_genus' in line:
            genus = pull_xml_record(line)
        elif "BinomialOrgName_species" in line:
            species = pull_xml_record(line)
        elif '<Dbtag_db>GeneID' in line:
            geneid_flag = 1
        elif geneid_flag == 1 and 'Object-id_id' in line:
            geneid = pull_xml_record(line)
        elif 'OrgName_lineage' in line:
            lineage = pull_xml_record(line)
            lineage = lineage.split(";")
            lineage = [x.strip() for x in lineage]
            lineage = ",".join(lineage)
        elif 'Gene-ref_locus-tag' in line:
            gene_tag = pull_xml_record(line)
        else:
            continue

