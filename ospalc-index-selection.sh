#Usage: bash ospalc-index-selection.sh forward-primer reverse-primer
#For example: bash ospalc-index-selection.sh GTGYCAGCMGCCGCGGTAA GGACTACNVGGGTWTCTAAT
#This script is used to rank the suitability of some variable components in the long primers, especially the 8nt indices in OSPALC primers. This script is based on short-blast. Depending on blast, selected forward and reverse indices are in forward/reverse.index.txt files, complete primers are in forward/reverse.primer.txt files, index filling rules for different sequencing platform are in forward/reverse.input.txt files.

a=AATGATACGGCGACCACCGAGATCTACAC #P5 sequence
b=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG #SP1 sequence
c=CAAGCAGAAGACGGCATACGAGAT #P7 sequence
d=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG #SP2 sequence
f=$1 #Forward target primer
r=$2 #Reverse target primer

p='{A,G,C,T}'
eval echo $p$p$p$p$p$p$p$p >index-8nt.txt #Generate all possible 8nt indices
sed -i "s/ /\n/g" index-8nt.txt
grep -E 'A|G|C|T' index-8nt.txt > indexAGCT.txt #Retain combinations that contain all four nucleotides
wc -l indexAGCT.txt > indexAGCT.number
cat indexAGCT.number
grep -v -E 'AAA|TTT|GGG|CCC' indexAGCT.txt > index-remove.txt #Retain indices without homopolymers ≥ 3nt
wc -l index-remove.txt > remove.number
cat remove.number
paste index-remove.txt index-remove.txt > primer-remove.fasta
sed -i "s/^/>/" primer-remove.fasta
sed -i "s/\t/\n/" primer-remove.fasta #Each index is presented in fasta format, titled with each index sequence
mkdir index-select
while read i; do grep -P -B 1 ^$i primer-remove.fasta > index-select/$i.fasta ; done < index-remove.txt #Split each sequence into a single file
cd index-select
while read i; do makeblastdb -in $i.fasta -dbtype nucl -out $i ; done < ../index-remove.txt 
mkdir ../index-select-blast/
while read i; do blastn -num_threads 8 -task blastn-short -word_size 4 -query $i.fasta -db $i -outfmt 7 -out ../index-select-blast/$i.txt ; done < ../index-remove.txt #Blast each index with itself in order to remove indices that are reverse complementary or similar to itself with -word_size=4
cd ../index-select-blast/
cat *.txt > all.txt
grep -B2 " 1 hits found" all.txt > all-grep.txt #Retain sequences only have 1 hit (itself)
sed -i "s/# Database: //" all-grep.txt
sed -i "s/#.*//" all-grep.txt
sed -i "s/-.*//" all-grep.txt
sed -i '/^[ ]*$/d' all-grep.txt #Remove the blast comments, retain only the sequence itself
wc -l all-grep.txt > ../remove-index.number
cat ../remove-index.number
mkdir ../preprimerf #Remove indices may cause self amplification in forward and reverse primers in directory preprimerf
cd ../preprimerf
cp ../index-select-blast/all-grep.txt ./
paste all-grep.txt all-grep.txt >preprimer.forward.fasta
sed -i "s/^/>/" preprimer.forward.fasta
sed -i "s/\t/\t$a/" preprimer.forward.fasta
sed -i "s/$/$b$f/" preprimer.forward.fasta
sed -i "s/\t/\n/" preprimer.forward.fasta #Generate fasta file by connecting the P5 sequence with the forward indices, the sequencing primer binding sequence and the forward target primer sequence, with the index sequence as the sequence title
while read i; do grep -P -B 1 ^$a$i preprimer.forward.fasta > $i.pre.forward.fasta ; done < all-grep.txt #Split each sequence into a single file
while read i; do makeblastdb -in $i.pre.forward.fasta -dbtype nucl -out $i ; done < all-grep.txt 
mkdir ../preprimerf_blast
while read i; do blastn -num_threads 8 -task blastn-short -word_size 4 -query $i.pre.forward.fasta -db $i -outfmt 7 -out ../preprimerf_blast/$i.pre.forward.txt ; done < all-grep.txt #Blast each forward primer with itself in order to remove primers that are reverse complementary or similar to itself with -word_size=4
cd ../preprimerf_blast
cat *.txt > all.pre.forward.txt
mkdir ../preprimerr
cd ../preprimerr
cp ../index-select-blast/all-grep.txt ./
paste all-grep.txt all-grep.txt >preprimer.reverse.fasta
sed -i "s/^/>/" preprimer.reverse.fasta
sed -i "s/\t/\t$c/" preprimer.reverse.fasta
sed -i "s/$/$d$r/" preprimer.reverse.fasta
sed -i "s/\t/\n/" preprimer.reverse.fasta #Generate fasta file by connecting the P7 sequence with the reverse indices, the sequencing primer binding sequence and the reverese target primer sequence, with the index sequence as the sequence title
while read i; do grep -P -B 1 ^$c$i preprimer.reverse.fasta > $i.pre.reverse.fasta ; done < all-grep.txt # Split each sequence into a single file
while read i; do makeblastdb -in $i.pre.reverse.fasta -dbtype nucl -out $i ; done < all-grep.txt 
mkdir ../preprimerr_blast 
while read i; do blastn -num_threads 8 -task blastn-short -word_size 4 -query $i.pre.reverse.fasta -db $i -outfmt 7 -out ../preprimerr_blast/$i.pre.reverse.txt ; done < all-grep.txt #Blast each reverse primer with itself in order to remove primers that are reverse complementary or similar to itself with -word_size=4
cd ../preprimerr_blast
cat *.txt > all.pre.reverse.txt
cd ..
grep hits preprimerf_blast/all.pre.forward.txt | sort -n -k2 | head -n 16 | uniq > preprimerf.hit
grep hits preprimerr_blast/all.pre.reverse.txt | sort -n -k2 | head -n 24 | uniq > preprimerr.hit #Get the least hit number in forward and reverse primers form the blast results above. In theory, at least 8 forward primers and 12 reverse primers are required, and 16 forward primers and 24 reverse primers are recommended here. Even if the 16 and 24 are set, the actual number per hit number will typically represent more than one sequence, and the output sequences will be much more than 16 and 24.
while read i ; do grep -B2 $i preprimerf_blast/all.pre.forward.txt >> preforward-primer.txt ; done < preprimerf.hit
sed -i "s/# Database: //" preforward-primer.txt
sed -i "s/#.*//" preforward-primer.txt
sed -i "s/-.*//" preforward-primer.txt
sed -i '/^[ ]*$/d' preforward-primer.txt #Remove the blast comments, only retain the selected forward indices 
wc -l preforward-primer.txt > preforward-primer.number
cat preforward-primer.number
while read i ; do grep -B2 $i preprimerr_blast/all.pre.reverse.txt >> prereverse-primer.txt ; done < preprimerr.hit
sed -i "s/# Database: //" prereverse-primer.txt
sed -i "s/#.*//" prereverse-primer.txt
sed -i "s/-.*//" prereverse-primer.txt
sed -i '/^[ ]*$/d' prereverse-primer.txt #Remove the blast comments, only retain the selected reverse indices 
wc -l prereverse-primer.txt > prereverse-primer.number
cat prereverse-primer.number
mkdir primerf
cd primerf
paste ../preforward-primer.txt ../preforward-primer.txt >primer.forward.fasta
sed -i "s/^/>/" primer.forward.fasta
sed -i "s/\t/\t$a/" primer.forward.fasta
sed -i "s/$/$b$f/" primer.forward.fasta
sed -i "s/\t/\n/" primer.forward.fasta #Generate fasta file by connecting the P5 sequence with the forward indices, the sequencing primer binding sequence and the reverese target primer sequence, with the index sequence as the sequence title
while read i; do grep -P -B 1 ^$a$i primer.forward.fasta > $i.forward.fasta ; done < ../preforward-primer.txt
makeblastdb -in primer.forward.fasta -dbtype nucl -out primerf.all
cd ..
mkdir primerr
cd primerr
paste ../prereverse-primer.txt ../prereverse-primer.txt > primer.reverse.fasta
sed -i "s/^/>/" primer.reverse.fasta
sed -i "s/\t/\t$c/" primer.reverse.fasta
sed -i "s/$/$d$r/" primer.reverse.fasta
sed -i "s/\t/\n/" primer.reverse.fasta #Generate fasta file by connecting the P7 sequence with the reverse indices, the sequencing primer binding sequence and the reverese target primer sequence, with the index sequence as the sequence title
while read i; do grep -P -B 1 ^$c$i primer.reverse.fasta > $i.reverse.fasta ; done < ../prereverse-primer.txt
makeblastdb -in primer.reverse.fasta -dbtype nucl -out primerr.all
while read i; do blastn -num_threads 8 -task blastn-short -word_size 4 -query $i.reverse.fasta -db ../primerf/primerf.all -outfmt 7 -out $i.reverse.txt ; done < ../prereverse-primer.txt #Blast each reverse primer with each forward primer in order to remove primers that are reverse complementary or similar with -word_size=4
mkdir ../primerr_blast
mv *.reverse.txt ../primerr_blast
cd ../primerr_blast
cat *.txt > all.reverse.txt
cd ../primerf
while read i; do blastn -num_threads 8 -task blastn-short -word_size 4 -query $i.forward.fasta -db ../primerr/primerr.all -outfmt 7 -out $i.forward.txt ; done < ../preforward-primer.txt #Blast each forward primer with each reverse primer in order to remove primers that are reverse complementary or similar with -word_size=4
mkdir ../primerf_blast
mv *.forward.txt ../primerf_blast
cd ../primerf_blast
cat *.txt > all.forward.txt
cd ..
grep hits primerf_blast/all.forward.txt | sort -n -k2 | head -n 16 | uniq > primerf.hit
grep hits primerr_blast/all.reverse.txt | sort -n -k2 | head -n 24 | uniq > primer.hit #Get the least hit number in forward and reverse primers form the blast results above. In theory, at least 8 forward primers and 12 reverse primers are required, and 16 forward primers and 24 reverse primers are recommended here. Even if the 16 and 24 are set, the actual number per hit number will typically represent more than one sequence, and the output sequences will be much more than 16 and 24.
while read i ; do grep -B3 $i primerf_blast/all.forward.txt >> forward.index.txt ; done < primerf.hit
sed -i "s/# Query: //" forward.index.txt
sed -i "s/#.*//" forward.index.txt
sed -i "s/-.*//" forward.index.txt
sed -i '/^[ ]*$/d' forward.index.txt #Remove the blast comments, only retain the selected forward indices 
wc -l forward.index.txt > forward.number 
while read i ; do grep -B3 $i primerr_blast/all.reverse.txt >> reverse.index.txt ; done < primer.hit
sed -i "s/# Query: //" reverse.index.txt
sed -i "s/#.*//" reverse.index.txt
sed -i "s/-.*//" reverse.index.txt
sed -i '/^[ ]*$/d' reverse.index.txt #Remove the blast comments, only retain the selected forward indices 
wc -l reverse.index.txt > reverse.number
comm forward.index.txt reverse.index.txt | cut -f 1 | sed '/^[ ]*$/d' > forward.comm-1.index.txt
cat forward.number
cat reverse.number
cat forward.comm-1.index.txt
paste forward.comm-1.index.txt forward.comm-1.index.txt > forward.primer.txt
sed -i "s/^/>/" forward.primer.txt
sed -i "s/\t/\t$a/" forward.primer.txt
sed -i "s/$/$b$f/" forward.primer.txt
sed -i "s/\t/\n/" forward.primer.txt #Generate selected complete forward primers in file forward.primer.txt
paste reverse.index.txt reverse.index.txt > reverse.primer.txt
sed -i "s/^/>/" reverse.primer.txt
sed -i "s/\t/\t$c/" reverse.primer.txt
sed -i "s/$/$d$r/" reverse.primer.txt
sed -i "s/\t/\n/" reverse.primer.txt #Generate selected complete reverse primers in file reverse.primer.txt while read i; do echo $i |tr a-z A-Z |tr ATCG TAGC |rev >> forward.comm-1.index-retr ; done < forward.comm-1.index.txt
paste forward.comm-1.index.txt forward.comm-1.index.txt forward.comm-1.index-retr | sed '1i Index\tNovaSeq 1.0、MiSeq、HiSeq 2000\/2500\tNovaSeq 1.5、MiniSeq、NextSeq、HiSeq 3000\/4000、HiSeq X' > forward.index.input.txt
paste reverse.index.txt reverse.index.txt reverse.index.txt | sed '1i Index\tNovaSeq 1.0、MiSeq、HiSeq 2000\/2500\tNovaSeq 1.5、MiniSeq、NextSeq、HiSeq 3000\/4000、HiSeq X' > reverse.index.input.txt #Generate indices filling methods for different Illumina sequencing platforms in forward.index.input.txt and reverse.index.input.txt
