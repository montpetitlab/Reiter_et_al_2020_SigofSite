# download gtf files
while read inline
do
wget $inline 
done < gtf_ftp.txt

# process gtf files
for infile in *gtf.gz
do
  j=$(basename $infile .gz)
  gunzip ${infile}
  tail -n +4 ${j} > ${j}_cut1
done

for infile in *cut1
do
  j=$(basename $infile .gtf_cut1)
  sed '$d' ${j}.gtf_cut1 >> all.gtf
done 

# convert MetFru2 from gff3 to gtf
# see R script metfru2_gff_to_gtf.R; removed first 3 lines by hand

# add genome name to augustus gtf gene names
# see R script augustus_gtf.R
# add in four extra gtfs from augustus and JGI annotations
for infile in GCA_002024125.1_ASM202412v1_genomic.gtf GCA_002894445.1_ASM289444v1_genomic.gtf GCA_002921095.1_ASM292109v1_genomic.gtf Metfru2_GeneCatalog_20170818.gtf
do
cat $infile >> all.gtf
done

gzip all.gtf

# download fasta files
while read inline
do
wget $inline 
done < fasta_ftp.txt


cat *fna.gz *fasta.gz > all.fna.gz

# Observations
# ANFW02000074.1 Metschnikowia fructicola 277 unitig_0, whole genome shotgun sequence	139
# APLS01000001.1 Hanseniaspora uvarum DSM 2768 Contig_1, whole genome shotgun sequence	116
# KL584974.1 Aureobasidium pullulans var. pullulans EXF-150 unplaced genomic scaffold scaffold_1, whole genome shotgun sequence	91
# MSDY01000045.1 Aureobasidium sp. FSWF8-4 FSW8-bin-4-scaffold_100, whole genome shotgun sequence	66
# ASM14353v4  Botrytis cinerea B05.10 (ascomycetes)	42
# MWSF01000001.1 Starmerella bacillaris strain FRI751 scaffold1, whole genome shotgun sequence	35
# LPNK01000115.1 Metschnikowia sp. AWRI3582 Contig_100, whole genome shotgun sequence	19
# CU928165.1 Lachancea thermotolerans CBS 6340 chromosome A complete sequence	12
# PEGC01000001.1 Cladosporium sp. SL-16 contig1, whole genome shotgun sequence	12
# LPNL01000001.1 Hanseniaspora opuntiae strain AWRI3578 Hanseniaspora_opuntiae_AWRI3578_scaffold1, whole genome shotgun sequence	8
# ALNQ01000001.1 Pichia kudriavzevii M12 contig1, whole genome shotgun sequence	4
# ASM129837v2  Saccharomyces sp. 'boulardii' (ascomycetes)	3
# JNDS01000001.1 Rhizopus stolonifer B9770 ctg7180000269960_1, whole genome shotgun sequence	3
# KB707673.1 Botryotinia fuckeliana BcDW1 unplaced genomic scaffold Scaffold_1, whole genome shotgun sequence	3
# JFAV02000001.1 Hanseniaspora vineae strain T02/19AF contig1, whole genome shotgun sequence	2
# LJJI01000001.1 Preussia sp. BSL10 Seq1, whole genome shotgun sequence	2
# LOQJ01000001.1 Saccharomyces pastorianus strain Hybrid yeast 1 sequence_1, whole genome shotgun sequence	2
# AABY01000832.1 Saccharomyces paradoxus NRRL Y-17217 contig_531, whole genome shotgun sequence	1
# FUGC01000105.1 Zygosaccharomyces bailii strain IST302 genome assembly, contig: ZBIST_scaffold106, whole genome shotgun sequence	1
# JFJX01000001.1 Tatumella sp. UCD-D_suzukii contig_1, whole genome shotgun sequence	1

# NO GTF
# organisms with JGI GTF files
https://genome.jgi.doe.gov/portal/Metfru2/download/Metfru2_AssemblyScaffolds.fasta.gz
https://genome.jgi.doe.gov/portal/Metfru2/download/Metfru2_GeneCatalog_20170818.gff3.gz
# GCA_000286515.1 Pichia kudriavzevii closest relative by sourmash k31 (88% similar):
GCA_003054445.1
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Pichia_kudriavzevii/latest_assembly_versions/GCA_003054445.1_ASM305444v1/GCA_003054445.1_ASM305444v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Pichia_kudriavzevii/latest_assembly_versions/GCA_003054445.1_ASM305444v1/GCA_003054445.1_ASM305444v1_genomic.gtf.gz
# Rhizopus stolonifer only/closest relative with GTF:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/325/415/GCA_003325415.1_Rstol_CA/GCA_003325415.1_Rstol_CA_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/325/415/GCA_003325415.1_Rstol_CA/GCA_003325415.1_Rstol_CA_genomic.gtf.gz
# Aureobasidium FSWF8 closet relative by sourmash k31 (90% similar):
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/917/305/GCA_004917305.1_ASM491730v1/GCA_004917305.1_ASM491730v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/917/305/GCA_004917305.1_ASM491730v1/GCA_004917305.1_ASM491730v1_genomic.gtf.gz
# MWSF01000001.1 Starmerella bacillaris strain FRI751 (model Saccharomyces cerevisiae [order match])
# annotated with WebAugustus 
# LPNK01000115.1 Metschnikowia sp. AWRI3582 no close relatives (model Saccharomyces cerevisiae [order match])
# annotated with WebAugustus
# PEGC01000001.1 Cladosporium sp. SL-16
# annotated with WebAugustus, model Aspergillus fumigatus [subphylum match] 

# Download links
# grape
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/745/GCA_000003745.2_12X/GCA_000003745.2_12X_genomic.fna.gz
# APLS01000001.1 Hanseniaspora uvarum DSM 2768
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/968/475/GCA_000968475.1_Huva2.0/GCA_000968475.1_Huva2.0_genomic.fna.gz
# ANFW02000074.1 Metschnikowia fructicola 277
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/317/355/GCA_000317355.2_ASM31735v2/GCA_000317355.2_ASM31735v2_genomic.fna.gz
# KL584974.1 Aureobasidium pullulans var. pullulans EXF-150 
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Aureobasidium_pullulans/latest_assembly_versions/GCA_000721785.1_Aureobasidium_pullulans_var._pullulans_EXF-150_assembly_version_1.0/GCA_000721785.1_Aureobasidium_pullulans_var._pullulans_EXF-150_assembly_version_1.0_genomic.fna.gz
# MSDY01000045.1 Aureobasidium sp. FSWF8-4 
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/914/275/GCA_001914275.1_FSWF8-bin-4/GCA_001914275.1_FSWF8-bin-4_genomic.fna.gz
# MWSF01000001.1 Starmerella bacillaris strain FRI751
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/024/125/GCA_002024125.1_ASM202412v1/GCA_002024125.1_ASM202412v1_genomic.fna.gz
# CU928165.1 Lachancea thermotolerans CBS 6340 c
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/142/805/GCF_000142805.1_ASM14280v1/GCF_000142805.1_ASM14280v1_genomic.fna.gz
# ALNQ01000001.1 Pichia kudriavzevii M12 
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/286/515/GCA_000286515.1_ASM28651v1/GCA_000286515.1_ASM28651v1_genomic.fna.gz
# KB707673.1 Botryotinia fuckeliana BcDW1
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Botrytis_cinerea/latest_assembly_versions/GCA_000349525.1_BcDW1/GCA_000349525.1_BcDW1_genomic.fna.gz
# ASM14353v4  Botrytis cinerea B05.10 
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/535/GCA_000143535.4_ASM14353v4/GCA_000143535.4_ASM14353v4_genomic.fna.gz
# LPNK01000115.1 Metschnikowia sp. AWRI3582 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/894/445/GCA_002894445.1_ASM289444v1/GCA_002894445.1_ASM289444v1_genomic.fna.gz
# PEGC01000001.1 Cladosporium sp. SL-16 
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/921/095/GCA_002921095.1_ASM292109v1/GCA_002921095.1_ASM292109v1_genomic.fna.gz
# LPNL01000001.1 Hanseniaspora opuntiae strain AWRI3578 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/749/795/GCA_001749795.1_ASM174979v1/GCA_001749795.1_ASM174979v1_genomic.fna.gz
# JNDS01000001.1 Rhizopus stolonifer B9770 
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/697/035/GCA_000697035.1_RhiStoB9770-1.0/GCA_000697035.1_RhiStoB9770-1.0_genomic.fna.gz
# JFAV02000001.1 Hanseniaspora vineae strain T02/19AF 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/585/475/GCA_000585475.3_FQ_Hvineae_v3/GCA_000585475.3_FQ_Hvineae_v3_genomic.fna.gz
# LJJI01000001.1 Preussia sp. BSL10 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/553/865/GCA_001553865.1_ASM155386v1/GCA_001553865.1_ASM155386v1_genomic.fna.gz

