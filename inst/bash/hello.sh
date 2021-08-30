###### Inputs:

echo -e "Working with RL: $RL"

###### Outputs:

outD="$2" # Specificed by USER: arg.2
BAMfile="$1"
fullname=$(basename $BAMfile)
RL=$(echo $fullname | cut -d . -f1)

dir4=$outD"/4.Alignment"
dir5=$outD"/5.Results"

mkdir -p $dir4 $dir5

###### Utilis:
samtools=$path_to_samtools
sambamba=$path_to_sambamba
macs2=$path_to_macs2
bamCoverage=$path_to_bamCoverage
bedtools=$path_to_bedtools
intersectBed=$path_to_intersectBed
mergeBed=$path_to_mergeBed

cp $BAMfile $dir4

###### Analyze:

# 5. STATS 给出BAM文件的比对结果
echo -e "\t5. Flag Stats"

$samtools flagstat $dir4"/"${RL}.bam  > $dir4"/"${RL}_flagStats
num_reads=$(grep QC $dir4"/"${RL}_flagStats | sed 's/+.*//g' )

# 6. Count MT reads
echo -e "\t6. Count MT reads"
$samtools index  $dir4"/"${RL}.bam
$samtools view -c $dir4"/"${RL}.bam chrM > $dir4"/"${RL}_mitoReads
mito_reads=$(cat $dir4"/"${RL}_mitoReads)

# 7. Proper Pairs
echo -e "\t7. Data is PE ---> $pairs\n\tFilter alignment, proper Pairs and stats"
$sambamba view -t 4 -f bam -F "proper_pair" $dir4"/"${RL}.bam -o $dir4"/"${RL}_pairs.bam

# 8. PCR dups
echo -e "\t8. Remove PCR dups and stats"
$sambamba markdup -p -t 4 -r --hash-table-size=1000000  $dir4"/"${RL}_pairs.bam $dir4"/"${RL}_pairs_dedup.bam

# 9. Remove blacklist
echo -e "\t9. Intercept BlackListed regions"
$samtools view -q 30 -L $blacklist -U $dir4"/"${RL}_pairs_dedup_filt.bam -b $dir4"/"${RL}_pairs_dedup.bam > $dir4"/"${RL}_pairs_dedup_crap

# 10. Remove MT reads
$samtools view -L $black_M -U $dir4"/"${RL}_pairs_dedup_filt_noMT.bam -b $dir4"/"${RL}_pairs_dedup_filt.bam > $dir4"/"${RL}_pairs_dedup_2crap
$samtools flagstat $dir4"/"${RL}_pairs_dedup_filt_noMT.bam > $dir4"/"${RL}_pairs_dedup_filt_noMT_flagStats
readsAfterFilter=$(grep QC $dir4"/"${RL}_pairs_dedup_filt_noMT_flagStats | sed 's/+.*//g' )

# Coverage --> may break the pipeline!
$samtools index $dir4"/"${RL}_pairs_dedup_filt_noMT.bam
$bamCoverage --binSize 50 -p 5 -b $dir4"/"${RL}_pairs_dedup_filt_noMT.bam -o $dir4"/"${RL}_pairs_dedup_filt_noMT_coverage.bigwig # maybe problems here, depending on size

# 11. RMs
echo -e "\t11. Remove Intermediate Files, some bam also"
rm $dir4"/"${RL}_*crap
rm $dir4"/"${RL}_pairs.bam*
rm $dir4"/"${RL}_pairs_dedup.bam*
rm $dir4"/"${RL}_pairs_dedup_filt.bam

# Just Keep original alignment + final filtered bam + index final + flagStats all + coverage.bw final
fBam=$dir4"/"${RL}_pairs_dedup_filt_noMT.bam

# 12. Reads in Promotors
echo -e "\t12. Counting reads in promoters"
$samtools view -c -L $promoters -b $fBam > $dir4/$RL"_ReadsInProm"
readsInPro=$(cat $dir4/$RL"_ReadsInProm")

# 13. MACS2 NFR peaks
echo -e "\t13. MACS2:::params for Tn5 insert sites, params ENCODE pipe for NFRs\n\tCalcule Reads in NFR peaks\n\tCalcule Reads in NFR peaks overlapping promotors"
$macs2 callpeak --nomodel --extsize 150 --shift -75 -t $fBam -f BAM -n $dir4"/"$RL"_NFR" --keep-dup all --gsize hs

cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | $samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInNFRpeaks"
readsInNFR=$(cat $dir4"/"$RL"_ReadsInNFRpeaks")

cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | $intersectBed -wa -a - -b $promoters | $samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInNFRpeaks_overProm"
readsInNFRoverPro=$(cat $dir4"/"$RL"_ReadsInNFRpeaks_overProm")

# 14. MACS BAMPE peaks
echo -e "\t14. MACS2:::params for PE fragments TLEN"
$macs2 callpeak -t $fBam -f BAMPE -n $dir4"/"$RL"_BAMPE" --keep-dup all --gsize hs
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | $samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInBAMPEpeaks"
readsInBAMPE=$(cat $dir4"/"$RL"_ReadsInBAMPEpeaks")

cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | $intersectBed -wa -a - -b $promoters | $samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBAMPEfrag_overProm"
readsInBAMPEoverPro=$(cat $dir4"/"$RL"_ReadsInBAMPEfrag_overProm")

$samtools view -f66 $fBam | cut -f 9 | sed 's/^-//' >  $dir4"/"$RL"_InsertSizesBAMPE"


# 15. Reads and Peaks In Background
echo -e "\t15. Reads in Background"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $bedtools subtract -wb -a $chrBed -b - | wc -l > $dir4"/"$RL"_PeaksInBackground

cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $bedtools subtract -a $chrBed -b - | $samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBackground"

readsInBack=$(cat $dir4"/"$RL"_ReadsInBackground")

echo "num_reads", "mito_reads", "readsAfterFilter", "readsInPro", "readsInNFR", "readsInNFRoverPro", "readsInBAMPE", "readsInBAMPEoverPro", "readsInBack" >  $dir4"/"$RL"_colnames.csv"
echo $num_reads, $mito_reads, $readsAfterFilter, $readsInPro, $readsInNFR, $readsInNFRoverPro, $readsInBAMPE, $readsInBAMPEoverPro, $readsInBack >  $dir4"/"$RL"_numreads.csv"

