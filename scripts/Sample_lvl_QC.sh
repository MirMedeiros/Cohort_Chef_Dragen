#!/bin/bash
#SBATCH --time=12:0:0
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --job-name=Sample_lvl_QC

################# SET UP VARIABLES #################

jointcalled_vcf=$1
genome_build=$4
genome_ref=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh${genome_build}/genome/Homo_sapiens.GRCh${genome_build}.fa
clinical_sexes=$2
output_dir=$3
bam_files=$5
script_dir=$PWD
P_DIR=$(echo "$PWD" | sed 's|/[^/]*$||')
cp ${P_DIR}/lib/all_hg${genome_build}.ref ${output_dir}
ref1KGPLINK=${P_DIR}/lib/1KG_hg${genome_build}_maf_hwe_geno

# make outdir and navigate to it:
mkdir ${output_dir}/Sample_QC
cd ${output_dir}/Sample_QC

########### Make some files for later:

##### Get Sample IDs
# load module
module load StdEnv/2023  gcc/12.3 bcftools/1.22

# extract sample IDs from VCF and save to file:
bcftools query -l $jointcalled_vcf > Sample_IDs.txt

#### Progress file set up:
touch SampleQC_progress.txt


###################################
#### CALCULATE VQSR ###############

### Calculate VQSR
###################
module purge
module load mugqic/java/openjdk-jdk-17.0.1 mugqic/GenomeAnalysisTK/4.6.0.0 mugqic/R_Bioconductor/4.1.0_3.13 && \
gatk --java-options "-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx24000M" \
  VariantRecalibrator  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh${genome_build}/genome/Homo_sapiens.GRCh${genome_build}.fa \
  --variant  ${jointcalled_vcf} \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.hapmap_3.3.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.1000G_omni2.5.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.1000G_phase1.snps.high_confidence.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.dbSNP142.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --truth-sensitivity-tranche 100.0 --truth-sensitivity-tranche 99.95 --truth-sensitivity-tranche 99.9 --truth-sensitivity-tranche 99.95 --truth-sensitivity-tranche 99.5 --truth-sensitivity-tranche 99.0 --truth-sensitivity-tranche 90.0 -mode SNP \
  --output allSamples.hc.snps.recal \
  --tranches-file allSamples.hc.snps.tranches && \
gatk --java-options "-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx24000M" \
  VariantRecalibrator  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh${genome_build}/genome/Homo_sapiens.GRCh${genome_build}.fa \
  --variant  ${jointcalled_vcf}  \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.Mills_and_1000G_gold_standard.indels.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh${genome_build}/annotations/Homo_sapiens.GRCh${genome_build}.dbSNP142.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --truth-sensitivity-tranche 100.0 --truth-sensitivity-tranche 99.9 --truth-sensitivity-tranche 99.4 --truth-sensitivity-tranche 99.0 --truth-sensitivity-tranche 90.0 -mode INDEL \
  --output allSamples.hc.indels.recal \
  --tranches-file allSamples.hc.indels.tranches

### Apply VQSR
##################

module load mugqic/java/openjdk-jdk-17.0.1 mugqic/GenomeAnalysisTK/4.6.0.0 && \
gatk --java-options "-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx24000M" \
  ApplyVQSR  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh${genome_build}/genome/Homo_sapiens.GRCh${genome_build}.fa \
  --variant ${jointcalled_vcf} \
  --truth-sensitivity-filter-level 99.95 --mode SNP \
  --tranches-file allSamples.hc.snps.tranches \
  --recal-file allSamples.hc.snps.recal \
  --output allSamples.hc.snps_raw_indels.vqsr.vcf.gz && \
gatk --java-options "-Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx24000M" \
  ApplyVQSR  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh${genome_build}/genome/Homo_sapiens.GRCh${genome_build}.fa \
  --variant allSamples.hc.snps_raw_indels.vqsr.vcf.gz \
  --truth-sensitivity-filter-level 99.4 --mode INDEL \
  --tranches-file allSamples.hc.indels.tranches \
  --recal-file allSamples.hc.indels.recal \
  --output allSamples.hc.vqsr.vcf.gz

rm  allSamples.hc.snps_raw_indels.vqsr.vcf.gz

####################################################
# Clean Up VCF (remove VQSR regions)

# load appropriate modules:
module load StdEnv/2023
module load gatk/4.4.0.0
module load python/3.13.2

gatk SelectVariants \
  -R $genome_ref \
  -V allSamples.hc.vqsr.vcf.gz \
  -O VQSR_Pass_vcf.gz \
  --java-options "-Xmx7g" \
  --tmp-dir ${SLURM_TMPDIR} \
  --exclude-filtered \
  --exclude-non-variants

## Update Progress:
echo "VSQR Filter Pass" SampleQC_progress.txt >> SampleQC_progress.txt

####################################################
############ Assemble QC Parameters: ###############


# Get GQ (Quality average)
# load module
module load StdEnv/2023  gcc/12.3 bcftools/1.22

### CALCULATE: #####################
####################################

# get GQ average:
bcftools query -f '[%SAMPLE\t%GQ\n]' "VQSR_Pass_vcf.gz" \
| awk -F'\t' '$2 != "." { sum[$1] += $2; n[$1]++ }
             END { for (s in sum) printf "%s\t%.6f\n", s, sum[s]/n[s] }' \
| sort -k1,1 > meanGQ.out

## Update Progress:
echo "Average Quality Calculated" SampleQC_progress.txt >> SampleQC_progress.txt


### CALCULATE CONTAMINATION AND CHIMERIC READS ###
##################################################

while read -r sample 
do

    # index cram file:
    module load StdEnv/2023 gcc/12.3 samtools/1.20
    samtools index ${sample}
    
    # extract name (accepts bams and crams):
    file_id=$(basename "$sample")
    ID=${file_id%.sorted.bam}
    ID=${ID%.sorted.cram}
    ID=${ID%.bam}
    ID=${ID%.cram}

    ##### get chimeric reads:
    module load picard/3.1.0

    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
              R=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh${genome_build}.fa  \
              I=${sample} \
              O=${ID}output_test.txt

    # extract chimeric reads value and save to temp variable:
    chimer=$(cut -f 1,29 ${ID}output_test.txt | grep -w "PAIR" | awk '{print $2}')

    ##### get contamination (source: https://github.com/Griffan/VerifyBamID):
    # get 1kG ref from github under resources
    module load mugqic/verifyBamID/2.0.1

    VerifyBamID \
      --UDPath ${VERIFYBAMID_HOME}/resource/1000g.phase3.10k.b${genome_build}.vcf.gz.dat.UD \
      --BedPath ${VERIFYBAMID_HOME}/resource/1000g.phase3.10k.b${genome_build}.vcf.gz.dat.bed \
      --MeanPath ${VERIFYBAMID_HOME}/resource/1000g.phase3.10k.b${genome_build}.vcf.gz.dat.mu \
      --Reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh${genome_build}.fa \
      --BamFile ${sample}

    # extract contamination value and save to temp variable:
    contam=$(awk '{print $7}' result.selfSM | tail -n 1)

    # save values
    echo $ID $chimer >> sample_chimeric.txt 
    echo $ID $contam >> sample_contamination.txt


done < ${bam_files}


### CALCULATE SAMPLE LVL DEPTH and MISSINGESS  ###
##################################################
module load StdEnv/2023  gcc/12.3
module load vcftools/0.1.16

VQSR_Pass_vcf.gz 

vcftools --gzvcf VQSR_Pass_vcf.gz --out depth_avg --depth
vcftools --gzvcf VQSR_Pass_vcf.gz --out missing_avg --missing-indv


##### PUT TOGETHER: ##################
######################################

## Update Progress:
echo "Chimeric Reads Extracted" SampleQC_progress.txt >> SampleQC_progress.txt

## Update Progress:
echo "Contamination Extracted" SampleQC_progress.txt >> SampleQC_progress.txt

########################################################
################ GENERATE PLINK FILES ##################

# load modules:
module load StdEnv/2023 plink/2.00-20231024-avx2

# start conversion - retaining only biallelic and autosomal SNPs:
plink2 --vcf VQSR_Pass_vcf.gz \
--double-id \
--max-alleles 2 \
--chr 1-22 \
--vcf-half-call missing \
--make-bed \
--memory 55000 \
--out VCF_PLINK

# Remain Vars to chr:pos:ref:alt
plink2 --bfile VCF_PLINK \
--memory 55000 \
--set-all-var-ids @:#:\$1:\$2 \
--new-id-max-allele-len 1000  \
--make-bed --out VarsRenamed

echo "." > missingSNPIDs.txt

plink2 --bfile VarsRenamed \
--memory 55000 \
--maf 0.05 \
--hwe 0.000001 midp \
--geno 0.02 \
--exclude missingSNPIDs.txt \
--make-bed \
--out Clean_VarsRenamed

rm VCF_PLINK*
rm missingSNPID.txt
rm VarsRenamed*

# Setup Fail Safe for if cohort is smaller than 30 samples to still run LD pruning:
# This will proceed with a sloppy higher noise pruning of the smaller cohort but be adaquet enough to look for general large differences between samples
SMALL_LD_FLAG=$( [ $(wc -l < Clean_VarsRenamed.fam) -lt 30 ] && echo "--bad-ld" )

plink2 --bfile Clean_VarsRenamed \
--indep-pairwise 50 5 0.5 ${SMALL_LD_FLAG} \
--memory 55000 \
--out Clean_VarsRenamed_pruning

plink2 --bfile Clean_VarsRenamed \
--exclude Clean_VarsRenamed_pruning.prune.out \
--memory 55000 \
--make-bed \
--out Pruned_PLINK

rm Clean_VarsRenamed_pruning*

# Update Progress
echo "Cohort PLINK Files Generated" SampleQC_progress.txt >> SampleQC_progress.txt

# run PCA:
plink2 \
--bfile Pruned_PLINK \
--pca 10 \
--out sample

# Update Progress
echo "Sample PCA Calculated" SampleQC_progress.txt >> SampleQC_progress.txt

##########################################################################
##
### Check Relatedness:
# second degree cut-off using the suggested threshold of 0.0884 in PLINK:

cleaned_PLINK=Clean_VarsRenamed

plink2 --bfile ${cleaned_PLINK} \
       --king-cutoff 0.0884 \
       --make-bed \
       --out unrelated_samples

# Related IDs are written here: unrelated_samples.king.cutoff.out.id

# Update Progress
echo "Relatedness Testing Done. Unrelated samples written to: unrelated_samples.king.cutoff.out.id" SampleQC_progress.txt >> SampleQC_progress.txt


##########################################################################
##
############# ANCESTRY TESTING #####################
# get list of overlapping SNPs
awk '{print $2}' ${cleaned_PLINK}.bim > sample_SNPs.txt
awk '{print $2}' ${ref1KGPLINK}.bim > ref_SNPs.txt

grep -f sample_SNPs.txt ref_SNPs.txt > overlapping_SNPs.txt

# get overlapping only plink sets:
plink2 --bfile ${cleaned_PLINK} \
--extract overlapping_SNPs.txt \
--make-bed \
--out overlapping_SNPs_samples

plink2 --bfile ${ref1KGPLINK} \
--extract overlapping_SNPs.txt \
--make-bed \
--out overlapping_SNPs_ref


# merge and QC together:
module load StdEnv/2020 plink/1.9b_6.21-x86_64

plink \
--bfile overlapping_SNPs_samples  \
--bmerge overlapping_SNPs_ref \
--make-bed \
--geno 0.001 \
--out merged_samples_ref

# clean
rm overlapping_SNPs*
rm sample_SNPs.txt
rm ref_SNPs.txt

# Prune Merged file:
module load StdEnv/2023 plink/2.00-20231024-avx2

plink2 --bfile merged_samples_ref \
--indep-pairwise 50 5 0.5 \
--memory 55000 \
--out samples_pruning

plink2 --bfile merged_samples_ref \
--exclude samples_pruning.prune.out \
--memory 55000 \
--make-bed \
--out pruned_merged

# Clean
rm samples_pruning*

# run PCA:
plink2 \
--bfile pruned_merged \
--pca 10 \
--out PCA_ancestry

# Update Progress
echo "Ancestry PCA Calculated" SampleQC_progress.txt >> SampleQC_progress.txt


##########################################################
###### CHECK SEX:
# keep sex chromosomes
plink2 --vcf VQSR_Pass_vcf.gz \
--double-id \
--max-alleles 2 \
--vcf-half-call missing \
--make-bed \
--memory 55000 \
--out joint_PLINK

# Remain Vars to chr:pos:ref:alt
plink2 --bfile joint_PLINK \
--memory 55000 \
--set-all-var-ids @:#:\$1:\$2 \
--new-id-max-allele-len 1000 missing \
--make-bed --out sex_chroms_incl_VarsRenamed

rm joint_PLINK*

# Do some standard QC of variants and exclude missing SNPIDs:
echo "." > missingSNPIDs.txt

plink2 --bfile sex_chroms_incl_VarsRenamed \
--memory 55000 \
--maf 0.05 \
--hwe 0.000001 midp \
--geno 0.02 \
--exclude missingSNPIDs.txt \
--make-bed \
--out Clean_cohort_sex_incl

# clean:
rm sex_chroms_incl_VarsRenamed*
rm missingSNPIDs.txt

# build up sexes files - 1 is male and 2 is female in PLINK:
awk '{print $1, $1, $2}' ${clinical_sexes} | sed 's/M/1/g' | sed 's/F/2/g' > sex_file.txt

# update sexes and check:
module load StdEnv/2020 plink/1.9b_6.21-x86_64

plink --bfile Clean_cohort_sex_incl \
      --update-sex sex_file.txt \
      --check-sex 0.4 0.6 \
      --out sex_check_output

# clean:
rm Clean_cohort_sex_incl*

# Identify samples which do not match with clinical sex:
grep "PROBLEM" sex_check_output.sexcheck > problem_samples.txt

# Update Progress
echo "Sex Check Complete. Problem samples writtent to: problem_samples" SampleQC_progress.txt >> SampleQC_progress.txt

###################################################################

module load  mugqic/R_Bioconductor/4.3.2_3.18

# Go into R and determine which Samples are bad and plot QC Metrics:
Rscript ${script_dir}/SampleQC.R "${output_dir}"

# Update Progress
echo "Identified Problematic Samples" SampleQC_progress.txt >> SampleQC_progress.txt

# Go into R and plot ancestries and cohort PCAs:
Rscript ${script_dir}/Ancestry_Plotting.R "${output_dir}" "${genome_build}"


# Update Progress
echo "Plotted Ancestry PCs" SampleQC_progress.txt >> SampleQC_progress.txt

###########################################

############### REMOVE FAILING SAMPLES ####################3
awk '{print $1}' unrelated_samples.king.cutoff.in.id > samples_unrelated.list 
grep -f samples_unrelated.list Samples_passing_QC.list > penultimate_passing_samples.txt

penultimate=$(wc -l penultimate_passing_samples.txt | awk '{print $1}')
echo "Unrelated Samples:" $penultimate >> Sample_QC_Passing_by_Step.txt

awk '{print $1}' problem_samples.txt > discordant_sex.txt  
grep -v -f discordant_sex.txt penultimate_passing_samples.txt > Final_list_of_passing_samples.txt

final_num=$(wc -l Final_list_of_passing_samples.txt | awk '{print $1}')
echo "Samples with Matching Sex:" $final_num >> Sample_QC_Passing_by_Step.txt

# list of samples which pass QC: Samples_passing_QC.list
sed '1d' Final_list_of_passing_samples.txt > passing_samples.txt

# filter:
bcftools view -S passing_samples.txt \
VQSR_Pass_vcf.gz \
 -o Post_Sample_QC_Post_VQSR.vcf.gz -O z

# Update Progress
echo "QC Failling Samples Removed" SampleQC_progress.txt >> SampleQC_progress.txt
echo "ALL SAMPLE QC STEPS COMPLETED! FINISHED! On to variant QC" SampleQC_progress.txt >> SampleQC_progress.txt

## CLEAN UP:
rm meanGQ.out sample_chimeric.txt sample_contamination.txt Pruned_PLINK.* Clean_VarsRenamed.* pruned_merged.* merged_samples_ref.* unrelated_samples.* penultimate_passing_samples.txt


echo "QC Step, Cut-off to Survive, Number of Remaining Samples" > temp
sed 's/Chimeric Pass: /Chimeric Reads, < 5%,/g' Sample_QC_Passing_by_Step.txt | sed 's/Contamination Pass: /Contamination, < 5%,/g' | sed 's/Missingness Pass: /Missingness, < Mean + 3 × SD,/g' | sed 's/Depth Pass: /Mean Depth, > Mean − 3 × SD,/g' | sed 's/Quality Pass: /Mean Quality, > Mean − 3 × SD,/g' | sed 's/Unrelated Samples: /Relatedness, Must have no relationship closer than 2nd degree,/g' | sed 's/Samples with Matching Sex: /Sex Check, Must not have sex mismatch,/g' > sqc.txt
cat temp sqc.txt | awk '!x[$0]++' > Sample_QC_Passing_by_Step.txt 


#### Move files to the output directory:
cp Ancestry_PCA_plot.png ${output_dir}
cp Cohort_PCs.png ${output_dir}
cp Cohort_screeplot.png ${output_dir}
cp discordant_sex.txt ${output_dir}
cp Final_list_of_passing_samples.txt ${output_dir}
cp sample.eigenvec ${output_dir}
cp Sample_QC_Passing_by_Step.txt ${output_dir}
cp Sample_QC_metrics.txt ${output_dir}
cp Samples_passing_QC.list ${output_dir}
cp PCA_ancestry.eigenval ${output_dir}
cp PCA_ancestry.eigenvec ${output_dir}
cp sample.eigenval ${output_dir}
cp qc_metrics_plot.png ${output_dir}
