#!/usr/bin/bash

#Creates a folder 
mkdir ~/Assignments/Assignment_1/Sequencing_data
#Copies all files ending with the file extension .fastq.gz from one directory to another
cp /media/Data_2/BIOL550/tmp/*.fastq.gz ~/Assignments/Assignment_1/Sequencing_data
#change directory
cd ~/Assignments/Assignment_1/Sequencing_data
#Decompresses all the .fastq.gz
gunzip -k *.fastq.gz
#Concatenates all the .fastq files into a single file
cat *.fastq > pacbio.fastq
#Changes the permissions of pacbio.fastq to 644
chmod 644 pacbio.fastq
echo 'PacBio data ready to be analyzed!'

