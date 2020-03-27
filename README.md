I'm testing the hypothesis that the Rhizipus genomes have different protein family domains that make them less susceptible to turbinmicin. 


First I pulled down all of the known sec14 genes from the ncbi to use as blast queries. Files and metadata: 
sequence/ncbi_sec14_genes.fa
sequence/ncbi_sec14_genes_taxonomy.xml
sequence/ncbi_sec14_genes.xml



```
## to pull down all of fungal genbank
ncbi-genome-download -s genbank -F cds-fasta -m genbank_metadata.txt -v fungi


## to pull down specific phyla from the genbank grabbed above 
~/scripts/tree_tools/fetch_refseq_taxonomy.py ~/sec14/db/mucoromycota genbank_metadata.txt False mucoromycota
~/scripts/tree_tools/fetch_refseq_taxonomy.py ~/sec14/db/zoopagomycota genbank_metadata.txt False zoopagomycota

## to grab the specific taxa used in the paper
ncbi-genome-download -s genbank -F cds-fasta -o sec14_genomes_genbank  -g "Candida auris,Candida albicans,Candida glabrata,Candida tropicalis,Aspergillus fumigatus,Fusarium,Scedosporium,Rhizopus" --flat-output -m sec14_genomes_metadata_genbank.txt -v fungi

```
I then copied all of the genomes of interest in to the 'all' directory: 'db/all'



Then I made blastdbs of all of these genomes, and blasted the known sec14 genes to them
```
cd db/all
ls | parallel makeblastdb -in {} -dbtype nucl -parse_seqids
mkdir tblastx

for i in `ls db/all/*fna`; 
	do base=`basename $i`; 
	tblastx -query sequence/ncbi_sec14_genes.fa -db $i -outfmt 6 -out tblastx/$base -evalue 0.001; 
	./scripts/parseBlast.R tblastx/$base tblastx/${base}.parsed; 
	done
	
```


Now subset out the sec14 CDS sequences from their main files.
```
for i in `ls tblastx/*parsed`; 
	do base=`basename $i`; 
	t=`cat $i`; 
	./scripts/subset_fasta.py db/all/${base%.parsed} $t >>output/sec14_all.fa; 
	done

```

And translate to amino acids
```
python3 ./scripts/translate_cds.py output/sec14_all.fa >output/sec14_all.aa
```

Run hmmer to identify protein coding domains in the sequences
```
hmmscan -o hmmer/sec14.pfam.out --domtblout hmmer/sec14.pfam.table -E 0.01 --domE 0.01 --cpu 10 /usr/local/db/pfam/Pfam-A.hmm output/sec14_all.aa 
awk '{print $4 "\t" $2 "\t" $3 "\t" $19 "\t" $20 "\t" $23 " " $24 " " $25}' hmmer/sec14.pfam.table  >hmmer/sec14.pfam.table.parsed
```


Make a tree of the sequences
```
~/scripts/tree_tools/make_tree.py output/sec14_all.fa 30 -a -cds
/usr/local/bin/raxmlHPC-PTHREADS-AVX -m GTRGAMMA -n all_sec14_cds_dialign  -s dialigntx/sec14_all.fa.phylip -f a -x 897543 -N autoMRE -p 345232 -T 40

python3 scripts/rename_newick.py RAxML_bipartitionsBranchLabels.all_sec14_cds_dialign output_renamed/sec14_all_id_map.txt
```

