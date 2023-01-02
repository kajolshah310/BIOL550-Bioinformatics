psswd
psswrd
passwrd
passwd
clear
ls -lh
ls -lha
mkdir Project_C
ls -lh
clear
ls -lh
cd Project_C/
clear
ls -lh
unzip 'Lists(1).zip' 
unzip scripts\(1\).zip 
clear
ls -lh
rm *.zip
clear
ls -lh
ls -lh Lists/
head Lists/SSRG.list 
clear
exit
clear
queryNCBI.pl 
which queryNCBI.pl 
clear
ls -lha
ls -lh /etc/inputrc 
nano /etc/inputrc 
nano .bashrc 
nano /etc/inputrc 
nano .bashrc 
source .bashrc
nano /etc/inputrc 
cd Project_C/
ls -llh
ls -lh Lists/
for each file List/*.csv; do echo $file; done
foreach file List/*.csv; do echo $file; done
for file List/*.csv; do echo $file; done
ls -lh
for file in List/*.csv; echo $file; done
for file in Lists/*.csv; echo $file; done
for file in Lists/*.csv; do echo $file; done
for file in Lists/*.csv \
for file in Lists/*.csv; do \
mkdir NCBI_DATA
cd NCBI_DATA/
for file in $PC/Lists/*.csv; do   queryNCBI.pl   -l $file   -fa   -gb   -p   -cds; done
for file in Lists/*.csv; do   queryNCBI.pl   -l $file   -fa   -gb   -p   -cds; done
for file in ../Lists/*.csv; do   queryNCBI.pl   -l $file   -fa   -gb   -p   -cds; done
clear
ls -lh
cd ..
clear
mkdir FastANI
cd FastANI/
ls ../NCBI_DATA/*.fasta > genome.list
head genome.list 
fastANI -t 6 -ql genome.list -rl genome.list -minFraction 0 -o fastANI.results.tsv
fastANI -t 6 --ql genome.list --rl genome.list --minFraction 0 -o fastANI.results.tsv
ls -lh
head fastANI.results.tsv 
ls -lh ../scripts/
chmod +x ../scripts/*
ls -lh ../scripts/
clear
head fastANI.results.tsv 
../scripts/remove_path.pl -in fastANI.results.tsv -out fastANI.clean.tsv -path ../NCBI_DATA/ 
head fastANI.results.tsv 
head fastANI.clean.tsv 
../scripts/fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
ls -lh
clear
MashR_plotter.pl -t heatmap -i ANI.tsv -if tsv -o ANI -f pdf -fs 8
MashR_plotter.pl
MashR_plotter.pl -l
ls -lh /
ls -lh /opt/
ls -lh /opt/mash-2.3/
clear
ls -lh
ls -lh fastANI.clean.tsv 
ls -lh ANI
nano ANI/Mash.R 
ANI/Mash.R 
ls -lh
nano ANI/Mash.R 
ls -lh
nano ../scripts/fastANItoMatrix.pl 
clear
cp ../scripts/fastANItoMatrix.pl ./
ls -lh
nano fastANI.clean.tsv 
ls -lh
cp fastANI.clean.tsv fastANI.clean.cp.tsv
nano fastANI.clean.cp.tsv 
nano fastANItoMatrix.pl 
./fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
clear
nano fastANItoMatrix.pl 
./fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
clear
./fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
clear
./fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
clear
nano fastANItoMatrix.pl 
./fastANItoMatrix.pl -i fastANI.clean.tsv -o ANI -f tsv
MashR_plotter.pl -t heatmap -i ANI.tsv -if tsv -o ANI -f pdf -fs 8
ls -lh
ls -lh ANI
ls -lh
rm fastANI.clean.cp.tsv 
clear
ls -lh
rm fastANItoMatrix.pl 
clear
ls -lh
ls -lh ANI
clear
ls -lh
ls -lh ../
fastANI -q ../NCBI_DATA/Clostridium_botulinum_Alaska_E43.fasta -r ../NCBI_DATA/Clostridium_botulinum_Mfbjulcb.fasta --visualize -o fastani.out
fastANI -q ../NCBI_DATA/Clostridium_botulinum_Alaska_E43.fasta -r ../NCBI_DATA/Clostridium_botulinum_Mfbjulcb8.fasta --visualize -o fastani.out
ls -lh
Rscript /opt/FastANI-1.33/scripts/visualize.R ../NCBI_DATA/Clostridium_botulinum_Alaska_E43.fasta ../NCBI_DATA/Clostridium_botulinum_Mfbjulcb8.fasta fastani.out.visual 
ls -lh
exit
mpirun -np Ray
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq SampleC_R2.nextera.fastq -o Ray_31
head SampleC_R1.nextera.fastq.gz 
clear
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq SampleC_R2.nextera.fastq -o Ray_31
ls -lh
ls 0lh Ray_31/
rm -r Ray_31/
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq SampleC_R2.nextera.fastq -o Ray_31
rm -r Ray_31/
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq SampleC_R2.nextera.fastq -o Ray_31
clear
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq.gz SampleC_R2.nextera.fastq.gz -o Ray_31
rm -r Ray_31/
mpirun -np 8 Ray -k 31 -p SampleC_R1.nextera.fastq.gz SampleC_R2.nextera.fastq.gz -o Ray_31
ls
ls Ray_31/
less Ray_31/LibraryStatistics.txt 
exitexit
exit
clear
ls -lh
cd Project_B/
ls -lh
fastqc
fastqc Alex_S4_L001_R*
ls -lh
fastqc --help
rm *.html
ls -lh
rm *.zip
clear
ls -lh
mv Alex_S4_L001_R1_001.fastq.gz SampleC_R1.fastq.gz
mv Alex_S4_L001_R2_001.fastq.gz SampleC_R2.fastq.gz
clear
ls -lh
mkdir BACKUP_READS
cp SampleC_R* BACKUP_READS/
clear
ls -lh
fastqc --help
fastqc *.fastq.gz -o FASTQC_PreFiltering
mkdir FASTQC_PreFiltering
fastqc *.fastq.gz -o FASTQC_PreFiltering
fastp- --help
fastp --help
clear
fastp --help
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapater_sequence_r2 CTGTCTCTTATACACATCT
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --threads 8
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8
ls -lh
mv fastp.* FASTQC_PreFiltering/
echo "fastqc *.fastq.gz -o FASTQC_PreFiltering" > genome_process.sh
echo "fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8" >> genome_process.sh 
clear
fastqc *nextera.fastq.gz -o FASTQC_PostFiltering
mkdir FASTQC_PostFiltering
fastqc *nextera.fastq.gz -o FASTQC_PostFiltering
ls
mv SampleC_R1.nextera.fastqc.gz SampleC_R1.nextera.fastq.gz
mv SampleC_R2.nextera.fastqc.gz SampleC_R2.nextera.fastq.gz 
fastqc *nextera.fastq.gz -o FASTQC_PostFiltering
echo "fastqc *nextera.fastq.gz -o FASTQC_PostFiltering" >> genome_process.sh 
mv Bacterium_C.pacbio.fastq.gz SampleC_pacbio.fastq.gz
fastqc SampleC_pacbio.fastq.gz -o FASTQC_PreFiltering/
ray
Ray
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100
echo "fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100" >> genome_process.sh 
fastqc *nextera.fastq.gz -o FASTQC_PostFiltering
clear
ls -lh
nano genome_process.sh 
mpirun
mpirun -np Ray -k 31 \
screen -S Ray
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100\
rm *nexter.fastq.gz
rm *nextera.fastq.gz
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100
rm *nextera*
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastqc.gz -O SampleC_R2.nextera.fastqc.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100
rm *nextera*
fastp -i SampleC_R1.fastq.gz -I SampleC_R2.fastq.gz -o SampleC_R1.nextera.fastq.gz -O SampleC_R2.nextera.fastq.gz --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --thread 8 -M 30 -r -l 100
fastqc *nextera.fastq.gz -o FASTQC_PostFiltering
screen -r
clear
ls -lh
screen -r
clear
exit
spades.py
cd Project_B/
ls -l
canu
exit
cd Project_
cd Project_B/
clear
ls -lh
dos2unix Genome_process.sh 
chmod +x Genome_process.sh 
screen -S assemble
exit
screen -r
exit
cd Project_B/
ls -lh
nano Genome_process.sh 
screen -r
spades.py
screen -r
exit
ls -lh 
cd Project_B/
ls -lh
ls -lh unicycler/
quast
quast.py
clear
ls -lh
mkdir Assemblies
cp Canu/Canu.contigs.fasta Assemblies/canu.fasta
clear
ls -lh
cp flye/assembly.fasta Assemblies/flye.fasta
clear
ls -lh
cp spades/contigs.fasta Assemblies/spades.fasta
clear
ls -lh
cp unicycler/assembly.fasta Assemblies/unicycler.fasta
quast.py
quast.py -o QUAST Assemblies/*.fasta
clear
ls -lh
ls -lh QUAST/
ls -lh /opt/
cd ..
ls -lh
git clone --recursive https://github.com/PombertLab/Misc.git
ls -lh Misc/
chmod +x Misc/*.pl
Misc/runTaxonomizedBLAST.pl 
Misc/runTaxonomizedBLAST.pl -t 4 -p blastn -a megablast -d nr -q Project_B/Assemblies/unicycler.fasta -e 1e-10 -c 1
Misc/runTaxonomizedBLAST.pl -t 4 -p blastn -a megablast -d NR -q Project_B/Assemblies/unicycler.fasta -e 1e-10 -c 1
less /etc/profile.d/bash.sh 
echo $BLASTDB
Misc/runTaxonomizedBLAST.pl -t 4 -p blastn -a megablast -d NR -q Project_B/Assemblies/unicycler.fasta -e 1e-10 -c 1
Misc/runTaxonomizedBLAST.pl -t 4 -p blastn -a megablast -d nr -q Project_B/Assemblies/unicycler.fasta -e 1e-10 -c 1
Misc/runTaxonomizedBLAST.pl -t 4 -q Project_B/Assemblies/unicycler.fasta -e 1e-10 -c 1
ls -lh /media/Data_1/NCBI/
ls -lh /media/Data_1/NCBI/NR/
ls -lh /media/Data_1/NCBI/NT
runTaxonomizedBLAST.pl    -t 6    -p blastn    -a megablast    -db ref_prok_rep_genomes    -q $PB/BLAST/*.fasta    -e 1e-30 \
Misc/runTaxonomizedBLAST.pl -t 6 -p blastn -a megablast -db ref_prok_rep_genomes -q Project_B/Assemblies/unicycler.fasta -e 1e-30 -c 1
ls -lh
ls -lh Project_B/
ls -lh
ls -lh Misc/
exit
ls -lh
ls -lh Project_B/
ls -lh Project_B/Assemblies/
less Project_B/Assemblies/unicycler.fasta.blastn.6 
Misc/runTaxonomizedBLAST.pl -t 6 -p blastn -a megablast -db ref_prok_rep_genomes -q Project_B/Assemblies/flye.fasta -e 1e-30 -c 1
less Project_B/Assemblies/flye.fasta.blastn.6 
grep -P "contig_2" Project_B/Assemblies/flye.fasta.blastn.6 
grep -P "contig_1" Project_B/Assemblies/flye.fasta.blastn.6 
grep -P "contig_2" Project_B/Assemblies/flye.fasta.blastn.6 
clear
exit
screen -S taxblast
exit
clera
clear
ls -lh
ls -lh Misc/
cd Project_B/
ls -lh
less Genome_process.sh 
history
fastp
fastp grep -P "-r"
grep -P "-r" fastp
fastp | grep -P "-r"
clear
fastp | grep -P "-r"
history
qless Genome_process.sh 
less Genome_process.sh 
which fastp
fastp
which fastqc
fastqc
spades.py
unicycler
unicycler --version
less Genome_process.sh 
history
quast.py --version
exit
cd Project_B/
ls -lh
fastp
ls -lh
fastp -i SampleC_pacbio.fastq.gz -l 500 -r -M 30
ls -lh
clear
blast
blastn
ls -lh
ls -lh Assemblies/
history
../Misc/parseTaxonomizedBLAST.pl 
../Misc/parseTaxonomizedBLAST.pl -b Assemblies/*.outfmt.6 \
cd Assemblies/
../../Misc/parseTaxonomizedBLAST.pl -b canu.fasta.blastn.6 -f canu.fasta \
head unicycler.fasta
head unicycler.fasta.blastn.6 
../../Misc/parseTaxonomizedBLAST.pl -b canu.fasta.blastn.6 -f canu.fasta -n "Vitreoscilla stercoraria" -out canu.parsed.fasta
ls -lh
../../Misc/parseTaxonomizedBLAST.pl -b canu.fasta.blastn.6 -f canu.fasta -n "Vitreoscilla stercoraria" -out canu.parsed.fasta -v
../../Misc/parseTaxonomizedBLAST.pl -b flye.fasta.blastn.6 -f flye.fasta -n "Vitreoscilla stercoraria" -out flye.parsed.fasta -v
../../Misc/parseTaxonomizedBLAST.pl -b spades.fasta.blastn.6 -f spades.fasta -n "Vitreoscilla stercoraria" -out spades.parsed.fasta -v
../../Misc/parseTaxonomizedBLAST.pl -b unicycler.fasta.blastn.6 -f unicycler.fasta -n "Vitreoscilla stercoraria" -out unicycler.parsed.fasta -v
clear
ls -lh
ls -lh ../
mv ../QUAST/ ./
clear
ls -lh
quast.py
quast.py *.parsed.fasta -o QUAST_PARSED
blastn --help
blastn -query unicycler.fasta -subject unicycler.fasta -o unicycler.compare
blastn -query unicycler.fasta -subject unicycler.fasta -out unicycler.compare
ls -lh
less unicycler.compare 
grep -P "Query=" unicycler.compare 
less unicycler.compare 
clear
ls -lh
prokka
cd ..
mv Assemblies/QUAST ./
mv Assemblies/QUAST_PARSED/ ./
prokka --addgenes --compliant --locustag BIO551 --increment 10 --genus Vitreoscilla --species stercoraria --gcode 11 --gram pos --cpus 12 Assemblies/unicycler.parsed.fasta 
ls -lh
ls -lh Assemblies/
ls PROKKA_03312022/
less PROKKA_03312022/PROKKA_03312022.txt 
python3 /opt/dfast_core-1.2.12/dfast --g Assemblies/unicycler.parsed.fasta --organism "Vitreoscilla stercoraria" --locus_tag_prefix BIOL551 --step 10 --out DFAST --cpu 12
/opt/dfast_core-1.2.12/dfast --g Assemblies/unicycler.parsed.fasta --organism "Vitreoscilla stercoraria" --locus_tag_prefix BIOL551 --step 10 --out DFAST --cpu 12
/opt/dfast_core-1.2.12/dfast --genome Assemblies/unicycler.parsed.fasta --organism "Vitreoscilla stercoraria" --locus_tag_prefix BIOL551 --step 10 --out DFAST --cpu 12
grep -c -P "^\s+(CDS)\s+" PROKKA_*/PROKKA_*gbk DFAST/gen*.gbk
ls -lh ../../Misc/
../Misc/runTaxonomizedBLAST.pl 
../../Misc/parseTaxonomizedBLAST.pl
../Misc/parseTaxonomizedBLAST.pl
exit
top
cd Project_B/
ls -lh
cd Assemblies/
ls -lh
grep flye.fasta.blastn.6 
less flye.fasta.blastn.6 
grep -P "contig_2" flye.fasta
grep -P "contig_2" flye.fasta.blastn.6 
grep -P "contig_2" canu.fasta.blastn.6 
grep -P -A20 "contig_2" flye.fasta
grep -P -A40 "contig_2"
grep -P -A40 "contig_2" flye.fasta
grep -P -A40 "contig_2" flye.fasta > flye.contig2_partial.fasta
blastn -query flye.contig2_partial.fasta -db nt -outfmt 6 -out flye.contig2_partial.fasta.blastn.6
blastn -query flye.contig2_partial.fasta -db nt -outfmt 6 -out flye.contig2_partial.fasta.blastn.6 -culling_limit 1
ls -lh
top
clear
ls -lh
nano flye.contig2_partial.fasta.blastn.6 
clear
ls -lh
nano ~/.bash_profile 
exit
ls -lh Project_B/
ls -lh Project_B/Ray_31/
less Project_B/Ray_31/CoverageDistributionAnalysis.txt 
less Project_B/Ray_31/CoverageDistribution.txt 
less Project_B/Ray_31/LibraryStatistics.txt 
clear
ls -lh
cd Project_B/
ls -lh
diff Ray_31/LibraryStatistics.txt Ray_alex/LibraryStatistics.txt 
clear
ls -lh
../Misc/keep_longest_reads.pl 
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -x -o Pacbio.metrics
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -m 1000 -o SampleC_pacbio.fastq.gz
git -C ../Misc/ pull
clear
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -m 1000 -o SampleC_pacbio.fastq.gz
ls -lh
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -x -o Pacbio.metrics
rm -r SampleC_pacbio.fastq.gz.log 
clear
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -x -o Pacbio.metrics
less SampleC_pacbio.fastq.gz 
mv SampleC_pacbio.fastq.gz SampleC_pacbio.fastq
pigz SampleC_pacbio.fastq 
ls -la
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -x -o Pacbio.metrics
ls -lh
ls -lh /media/
ls -lh /media/Data_3/
ls -lh /media/Data_3/jpombert/
ls -lh
rm -r SampleC_pacbio.fastq.gz 
ls -lh
mv Bacterium_C.pacbio.fastq.gz SampleC_pacbio.fastq.gz
ls -lh
rm -r fastp.html 
rm -r fastp.json 
rm -r metrics.log 
../Misc/keep_longest_reads.pl -i SampleC_pacbio.fastq.gz -m 1000 -o SampleC_pacbio_1k.fastq
ssh julian@216.47.151.157
clear
cd ..
ls -lh
cd Project_A/
ls -lha
ls -lh
cd ..
git clone --recursive https://github.com/PombertLab/ID16S.git
cd Project_A/
ls -lh
tar -xvf 16S_barcoded_S22.tar.gz 
clear
ls -lh
ls -lh 16S_
ls -lh 16S_Basecall_2022_03_28/
ls -lh 16S_Basecall_2022_03_28/pass/
cat 16S_Basecall_2022_03_28/pass/barcode11/*.fastq > barcode11.fastq
ls -lh
../ID16S/run_ID16S.pl
../ID16S/setup_ID16S.pl 
../ID16S/download_ID16S_dbs.pl 
../ID16S/setup_ID16S.pl 
../ID16S/setup_ID16S.pl -c ID16S.config -d ../ID16S/ -i ../ID16S/
source ID16S.config 
../ID16S/download_ID16S_dbs.pl 
../ID16S/download_ID16S_dbs.pl -d
../ID16S/run_ID16S.pl 
../ID16S/run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq -hd 15 -m 1000
screen -S blast
clear
./Genome_process.sh 
clear
nano Genome_process.sh 
./Genome_process.sh 
unicycler --help
unicycler -1 SampleC_R1.nextera.fastq.gz -2 SampleC_R2.nextera.fastq.gz -l SampleC_pacbio.fastq.gz -t 8 --pilon_path /opt/pilon/pilon-1.24.jar --spades_path /opt/SPAdes-3.13.0-Linux/
unicycler -1 SampleC_R1.nextera.fastq.gz -2 SampleC_R2.nextera.fastq.gz -l SampleC_pacbio.fastq.gz -t 8 --pilon_path /opt/pilon/pilon-1.24.jar --spades_path /opt/SPAdes-3.13.0-Linux/ -o unicycler
unicycler -1 SampleC_R1.nextera.fastq.gz -2 SampleC_R2.nextera.fastq.gz -l SampleC_pacbio.fastq.gz -t 8 --pilon_path /opt/pilon/pilon-1.24.jar --spades_path /opt/SPAdes-3.13.0-Linux/bin/ -o unicycler
unicycler -1 SampleC_R1.nextera.fastq.gz -2 SampleC_R2.nextera.fastq.gz -l SampleC_pacbio.fastq.gz -t 8 --pilon_path /opt/pilon/pilon-1.24.jar --spades_path /opt/SPAdes-3.13.0-Linux/bin/spades.py -o unicycler
clear
./Genome_process.sh 
nano Genome_process.sh 
./Genome_process.sh 
spades.py
spades.py | grep "--pe-1"
spades.py | grep -P "--pe-1"
grep -P "--pe-1" `spades.py`
spades.py
clear
nano Genome_process.sh 
./Genome_process.sh 
ls -lh
rm -r spades/
nano Genome_process.sh 
./Genome_process.sh 
nano Genome_process.sh 
./Genome_process.sh 
nano Genome_process.sh 
./Genome_process.sh 
nano Genome_process.sh 
./Genome_process.sh 
ls -lh
rm -r Canu/
nano Genome_process.sh 
./Genome_process.sh 
exit
Misc/runTaxonomizedBLAST.pl -t 6 -p blastn -a megablast -db ref_prok_rep_genomes -q Project_B/Assemblies/*.fasta -e 1e-30 -c 1
exit
ls -lh
../ID16S/run_ID16S.pl 
rm -r ID16S_RESULTS/
../ID16S/run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq 
ls -lh
ls -lh ID16S_RESULTS/
ls -lh ID16S_RESULTS/Normalized/
less ID16S_RESULTS/Normalized/barcode11_fasta_megablast_genus_Normalized_Microbiome_Composition.tsv 
clear
exit
ls -lh
screen -r
screen -r assemble 
screen -r taxblast 
screen -r blast 
cd Project_A
clear
ls -lh
less ID16S_RESULTS/NonNormalized/
less ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.class 
less ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species 
../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
ls -lh
ls -lh Normalized/
less Normalized/barcode11_fasta_megablast_species_Normalized_Microbiome_Composition.tsv 
sftp julian@216.47.151.157
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample bc04_16S.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
ls -lh Normalized/
less Normalized/bc04_16S_fasta_megablast_species_Normalized_Microbiome_Composition.tsv 
exit
cd Project_A
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
clear
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o Normalized
nano ../ID16S/Normalization_scripts/get_organism_statistics.pl 
exit
ls -lh /opt/
git -C /opt/ID16S/ pull
sudo git -C /opt/ID16S/ pull
ssh julian@216.47.151.148
clear
cd Project_A/
ls -lh
ls ID16S_RESULTS/
rm -r Normalized/
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.* -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
git -C ../ID16S/ pull
../ID16S/Normalization_scripts/get_organism_statistics.pl --sample ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.* -d ../ID16S/Normalization_DB/ -n ../ID16S/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
../ID16S/Normalization_scripts/get_organism_statistics.pl 
../ID16S/run_ID16S.pl 
../ID16S/run_ID16S.pl -o ID16S_RESULTS/ -fq barcode11.fastq
screen -S blast
clear
exit
ls -lh
history
clear
cd Project_B/
ls -lh
nano Genome_process.sh 
history
clear
history | grep "quast"
history | grep "ray"
history | grep "Ray"
which ray
which Ray
ls -lh /opt/Ray/
ls -lh /opt/Ray/README.md 
less /opt/Ray/README.md 
less Ray_alex/RayVersion.txt 
fastqc --help
fastqc --version
canu --version
clear
ls -lh
history
cd Project_
cd Project_B
ls -lh
cd Assemblies/
ls -lh
nano split_fasta.pl
chmod +x split_fasta.pl 
./split_fasta.pl flye.fasta
ls -lh
rm -r Split_Fasta/
./split_fasta.pl
./split_fasta.pl -f flye.fasta -o FLYE
ls -lh
ls -lh FLYE/
blastn -query FLYE/contig_1.fasta -subject unicycler.fasta -out uni_v_flye.blast -outfmt 6
ls -lh
less uni_v_flye.blast 
rm -r uni_v_flye.blast 
blastn -query FLYE/contig_1.fasta -subject unicycler.fasta -out uni_v_flye.blast -outfmt 0 -culling_limit 1
less uni_v_flye.blast 
cd ..
ls -lh
ls -lh PROKKA_03312022/
less PROKKA_03312022/PROKKA_03312022.txt 
less DFAST/
 DFAST/statistics.txt 
ls -lh DFAST/
head DFAST/cds.fna 
grep -pc ">" DFAST/cds.fna 
grep -Pc ">" DFAST/cds.fna 
grep -Pc "hypothetical" DFAST/cds.fna 
grep -Pc "hypothetical" PROKKA_03312022/PROKKA_03312022.fna 
grep -Pc "hypothetical" PROKKA_03312022/PROKKA_03312022.tsv
grep -Pc "hypothetical" DFAST/cds.fna 
grep -Pc "hypothetical" DFAST/protein.faa 
cd Project_B/
ls -lh
ls -lh Assemblies/
less Assemblies/unicycler.compare 
exit
../ID16S/run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq \
rm -r ID16S_RESULTS/
run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq \
run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq
exit
ls -lh
cd Project_A/
ls -lh
screen -r
rm -r ID16S_RESULTS/
ls -lh /opt/ID16S/ID16S_DB/Normalization_DB/
exit
clear
cd Project_A
screen -S ID16S
ls -lh
rm bc04_16S.fasta.megablast.species 
ls -lh ../
rm -r ../ID16S/
rm -rf ../ID16S/
clear
ls -lh
ls -lh ID16S_RESULTS/
ls -lh ID16S_RESULTS/NonNormalized/
ls -lh ID16S_RESULTS/Normalized/
ls -lh ID16S_RESULTS/FASTA/
exit
screen -r
exit
ls -lh
screen -r
exit
screen -r
exit
kraken2 --threads 4 --db 16S_RDP_k2db --report barcode_11_report.txt --report-minimizer-data --output barcode_11_output.txt barcode11.fastq 
find
clear
ls -lh
ls -lh /
ls -lh /media/
ls -lh /media/Data_1
which kraken2
ls -lh /opt/kraken2/
ls -lh /opt/kraken2/databases/
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_11_report.txt --report-minimizer-data --output barcode_11_output.txt barcode11.fastq 
ls -lh
less barcode_11_output.txt 
less barcode_11_report.txt 
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_11_report.txt --report-minimizer-data --output barcode_11_output.txt barcode11.fastq --use-names
less barcode_11_report.txt 
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_11_report.txt --report-minimizer-data --output barcode_11_output.txt barcode11.fastq
less barcode_11_report.txt 
exit
top
clear
cd Project_A
ls -lh
ls -lh barcode01/
ls -lh ID16S_BC1/
ls -lh ID16S_BC1/Normalized/
diff ID16S_BC1/Normalized/barcode1_fasta_megablast_order_Normalized_Microbiome_Composition.tsv ID16S_RESULTS/Normalized/barcode11_fasta_megablast_order_Normalized_Microbiome_Composition.tsv 
ls -lh
ls -lh Data_from_previous_year/
screen -r
screen -S KRAKEN
less ID16S_RESULTS/NonNormalized/ba
head ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.genus 
less ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.genus 
less ID16S_RESULTS/Normalized/barcode11_fasta_megablast_genus_Normalized_Microbiome_Composition.tsv 
ls -lh
ls -lh barcode_11_output.txt 
less barcode_11_output.txt 
less barcode_11_report.txt
exit
cd Project_
cd Project_A/
ls -lh
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_11_report.txt --output barcode_11_output.txt barcode11.fastq
less barcode_11_report.txt
dos2unix parse_kraken2_results.pl 
chmod +x parse_kraken2_results.pl 
exit
cd Project_A/
which kraken2
kraken2 --version
exit
top
clera
clear
ls -lh
cd Project_A/
ls -lh
ls -lh ID16S_BC1/
ls -lh ID16S_BC1/No
ls -lh ID16S_BC1/Normalized/
ls -lh ID16S_BC1/Normalized/Figure_files/
less /opt/ID16S/Normalization_scripts/create_figure.py 
exit
kraken2
cd Project_A/
exit
cd Project_
cd Project_A
kraken2
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_11_report.txt --output barcode_11_output.txt barcode11.fastq --use-names
less barcode_11_report.txt 
exit
cd Project_A/
ls -lh
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_1_report.txt --output barcode_1_output.txt barcode1.fastq --use-names
ls -lh
sftp julian@216.47.151.157
exit
cd Project_A
ls -lh
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_07_report.txt --output barcode_07_output.txt barcode07.fastq
sftp julian@216.47.151.157
exit
cd Project_A/
ls -lh
exit
cd Project_A
ls -lh
ls -lh Data_from_previous_year/
tar -xvf Data_from_previous_year/barcode.08.tar.gz 
ls -lh
screen -r
exit
ls -lh
screen -r
exit
screen -r
exit
run_ID16S.pl -o ID16S_RESULTS -fq barcode11.fastq
ls -lh
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl 
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl \
less /opt/ID16S/Normalization_scripts/get_organism_statistics.pl 
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.species -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.genus -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.family -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.class -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.family -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESULTS/Normalized/
/opt/ID16S/Normalization_scripts/get_organism_statistics.pl -s ID16S_RESULTS/NonNormalized/barcode11.fasta.megablast.order -d /opt/ID16S/ID16S_DB/Normalization_DB/ -n /opt/ID16S/ID16S_DB/TaxDump/names.dmp -o ID16S_RESU/Normalized/
ls -lh ID16S_RESULTS/Normalized/
less ID16S_RESULTS/Normalized/barcode11_fasta_megablast_class_Normalized_Microbiome_Composition.tsv 
less ID16S_RESULTS/Normalized/barcode11_fasta_megablast_species_Normalized_Microbiome_Composition.tsv 
clear
less ID16S_RESULTS/Normalized/barcode11_fasta_megablast_species_Normalized_Microbiome_Composition.tsv 
/opt/ID16S/run_ID16S.pl 
/opt/ID16S/run_ID16S.pl \
ls -lh Data_from_previous_year/
tar -xvf Data_from_previous_year/barcode.01.tar.gz 
ls -lh
cat Data_from_previous_year/ > barcode1.fastq
/opt/ID16S/run_ID16S.pl 
/opt/ID16S/run_ID16S.pl -o ID16S_BC1 -fq barcode1.fastq \
tar -xvf Data_from_previous_year/barcode.07.tar.gz 
cat barcode07/*.fastq > barcode07.fastq
/opt/ID16S/run_ID16S.pl 
/opt/ID16S/run_ID16S.pl -o ID16S_bc7 -fq barcode07.fastq -t 4
ls -lh
ls barcode08/
cat barcode08/*.fastq > barcode_08.fastq
/opt/ID16S/run_ID16S.pl -o ID16S_bc8 -fq barcode_08.fastq -t 4
exit
screen -r
exit
cd Project_A
kraken2 --threads 4 --db /opt/kraken2/databases/16S_RDP_k2db/ --report barcode_08_report.txt --output barcode_08_output.txt barcode_08.fastq 
exit
ls -lh
cd Project_C/
ls -lh
ls -lh FastANI/
history | grep 'ani'
ls -lh
ls -lh Lists/
ls -lh NCBI_DATA/
clear
ls -lh
ls -lh Lists/
less Lists/botulinum.csv 
clear
ls -lh
ls -lh FastANI/
ls -lh
ls -lh FastANI/
which fastani
which FastANI
clear
exit
cd Project_A/
ls -lh
cd ../Project_C/
ls -lh
clear
ls -lh FastANI/
ls -lh FastANI/ANI
less FastANI/ANI.tsv 
ls -lh FastANI/fastANI.results.tsv 
less FastANI/fastANI.results.tsv 
quast.py
which quast.py
python3 /opt/quast-5.0.2/quast.py 
ls -lh /opt
exit
