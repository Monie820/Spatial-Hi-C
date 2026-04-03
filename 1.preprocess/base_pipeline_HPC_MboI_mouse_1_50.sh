# source /share/home/zhanglin/.bashrc
. /date/scripts/Pipelines/temp_bashrc

### Software ###
fastqc=/share/home/share/software/FastQC/bin/fastqc
java_jre=/share/home/share/software/java/jre1.8.0_261/bin/java
cutadapt=/share/home/share/software/Python-3.9.2/bin/cutadapt
umi_tools=/share/home/zhanglin/softwares/miniconda3/bin/umi_tools
bwa=/share/home/share/software/bwa-0.7.17/bin/bwa
bedtools=/share/home/share/software/bedtools2/bin/bedtools
python_exec=/share/home/zhanglin/softwares/miniconda3/bin/python
r_exec=/share/home/zhanglin/softwares/miniconda3/bin/Rscript
hicFindRestSite=/share/home/zhanglin/softwares/miniconda3/envs/hicexplorer/bin/hicFindRestSite
hicBuildMatrix=/share/home/zhanglin/softwares/miniconda3/envs/hicexplorer/bin/hicBuildMatrix

ppn=16
core=4
sample=$2    # sample=Spatial-E13-20-20220127  sample=Spatial-E13-10-20220817   Spatial-Cere-38-20220926
work_path=$1/$sample    # work_path=/media/maoni/data/CZP/spatial_hic/Spatial-E13-20-20220127   work_path=/media/maoni/data/CZP/spatial_hic/Spatial-E13-10-20220817
script_path=$1/scripts_50        # script_path=/media/maoni/data/CZP/spatial_hic/scripts_50

data_path=$work_path/data
barcode_ref=$work_path/barcode_ref
check_linker_path=$work_path/0_check_linker
trimmed_path=$work_path/1_cutadapt
debarcode_path=$work_path/2_debarcode
hicMatrix_path=$work_path/3_hicMatrix
bam_path=$work_path/4_bam
bedpe_path=$work_path/5_bedpe

fa_path=/date/lvjunjie/Reference/mouse/GRCm38_mm10/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
cutsite=/date/guomaoni/Reference/HiC_files/GRCm38_MboI_rest_site_positions.bed
index=/date/lvjunjie/Reference/mouse/GRCm38_mm10/ensembl/bwa_index/GRCm38


if [ -d $work_path ] && [ ! -z $sample ]; then
	### fastqc ###  cat E13-9-20220701/E13-9-RNA-20220621_R1.fq.gz E13-9-20220714/E13-9-RNA-20220621_R1.fq.gz > E13-9-RNA-20220621_R1.fq.gz
	mkdir $work_path/fastqc
	$fastqc -t $ppn ${data_path}/${sample}_R1.fq.gz ${data_path}/${sample}_R2.fq.gz -o $work_path/fastqc
	
	
    ### 0.check linker (QC step1) ###
	mkdir $check_linker_path
	echo "Starting to process Hi-C data:" >> ${check_linker_path}/check_linker.log
	echo "Raw reads in total:" >> ${check_linker_path}/check_linker.log
	zcat ${data_path}/${sample}_R1.fq.gz | awk '{if (NR%4==2){print $0}}' - | wc -l >> ${check_linker_path}/check_linker.log  # Total reads: 485962701
	# check >=2 linkers
	echo "Multiple barcode1:" >> ${check_linker_path}/check_linker.log
	zcat ${data_path}/${sample}_R1.fq.gz | grep "CATCGGCGTACGACT[ATCG]\{22,24\}CATCGGCGTACGACT" | wc -l >> ${check_linker_path}/check_linker.log  #  Multiple barcode1: 15978866
	echo "Multiple barcode2:" >> ${check_linker_path}/check_linker.log
	zcat ${data_path}/${sample}_R1.fq.gz | grep "GTGGCCGATGTTTCG[ATCG]\{22,24\}GTGGCCGATGTTTCG" | wc -l >> ${check_linker_path}/check_linker.log
	echo "ME barcode1 mix:" >> ${check_linker_path}/check_linker.log
    zcat ${data_path}/${sample}_R1.fq.gz | grep "CATCGGCGTACGACT[ATCG]\{27,29\}CATCGGCGTACGACT" | wc -l >> ${check_linker_path}/check_linker.log
    zcat ${data_path}/${sample}_R1.fq.gz | awk '{if (NR%4==2){print $0}}' - | grep "CATCGGCGTACGACT[ATCG]\{22,24\}CATCGGCGTACGACT" > ${check_linker_path}/1_2ormore_barcodeA.txt
          
          
	### 1.trimmomatic ###
	mkdir $trimmed_path
	$cutadapt -j 6 \
    	         -a 'CTGTCTCTTATACACATCT' \
    	         -m 37 \
	             --trim-n \
    	         --pair-filter=first \
    	         -o ${trimmed_path}/${sample}.cutadapt.R2.fq.gz \
	             -p ${trimmed_path}/${sample}.cutadapt.R1.fq.gz \
	             ${data_path}/${sample}_R2.fq.gz ${data_path}/${sample}_R1.fq.gz
	
	### fastqc ###
	$fastqc -t $ppn ${trimmed_path}/${sample}.cutadapt.R1.fq.gz ${trimmed_path}/${sample}.cutadapt.R2.fq.gz -o $work_path/fastqc
	
	
	### 2.debarcode ###
	# Insertion, indicated by "i" (extra)
	# Deletion, indicated by "d" (missing)
	# Substitution, indicated by "s" (mismatch)
	# The above 3 errors, indicated by "e" (mismatch)
	mkdir $debarcode_path
	$umi_tools extract \
	--extract-method=regex \
	--bc-pattern="^(?P<cell_1>.{8})(?P<discard_1>GTGGCCGATGTTTCGCATCGGCGTACGACT){e<=1}(?P<umi_1>.{8})(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTGAGATGTGTATAAGAGACAG){s<=2,i<=1,d<=1}.*" \
	-I ${trimmed_path}/${sample}.cutadapt.R1.fq.gz -S ${debarcode_path}/${sample}.extract.R1.fq.gz \
	--read2-in=${trimmed_path}/${sample}.cutadapt.R2.fq.gz --read2-out=${debarcode_path}/${sample}.extract.R2.fq.gz \
	-L ${debarcode_path}/extract.log
          
    ### fastqc ###
	$fastqc -t $ppn ${debarcode_path}/${sample}.extract.R1.fq.gz ${debarcode_path}/${sample}.extract.R2.fq.gz -o $work_path/fastqc
          
	### 3. extract barcode to a txt ###
	zcat ${debarcode_path}/${sample}.extract.R1.fq.gz | awk 'NR%4==1{print $0}' | sed 's/^.*_\([A-Z]\{8\}\)_\([A-Z]\{8\}\) .*$/\1_\2/g' - > ${debarcode_path}/barcodeB_A.read1.txt
             
    mkdir $barcode_ref
    $r_exec $script_path/barcode_file_pre.R \
	-p $barcode_ref \
	-a $barcode_ref/barcodes_A.txt \
	-b $barcode_ref/barcodes_B.txt
	
    $r_exec $script_path/9_check_barcode.R \
	-p $check_linker_path \
	-n $debarcode_path/barcodeB_A.read1.txt \
	-b $barcode_ref/2_combine_barcode.round2round1_index1_index2.txt
	
	# "total barcode num"
          # 454584739   (454584739/485962701==ligation)
 	# "our barcode num"
          # 444126952   (444126952/454584739==barcode)

	$r_exec $script_path/check_multiple_barcodeA.R \
	-q $ppn \
	-d $check_linker_path \
	-i $check_linker_path/1_2ormore_barcodeA.txt \
	-b $barcode_ref/barcodes_A.txt \
	-o $barcode_ref/1_barcode_add_order.n1_n2.txt
	
          # calculate the 2,3 colums for the cross channel number: 2247739 + 7664 = 2255403 (2255403/15978866) 
          # unique.ls
          #          1        2        3 
          #    13723463  2247739     7664
	
	### 4. visualize pixel seq-stat ###
	$python_exec $script_path/check_barcode_visualize.py \
	--bc=$debarcode_path/barcodeB_A.read1.txt \
	--refA=$barcode_ref/barcodes_A.txt \
	--refB=$barcode_ref/barcodes_B.txt \
	--fig=$debarcode_path/pixel_rawcount_heatmap.png \
	--map_stat=$debarcode_path/raw_stat.csv \
	--unmap_stat=$debarcode_path/bad_pixel_stat.csv

	$r_exec $script_path/draw_svg.R \
	-f $debarcode_path/raw_stat.csv \
	-p $work_path
	
	### 5.HiCExplorer ###
	mkdir $bam_path
	$bwa mem -A 1 -B 4 -E 50 -L 0 -t $ppn ${index} ${debarcode_path}/${sample}.extract.R1.fq.gz | samtools view -Shb - > ${bam_path}/${sample}_R1.bam
	$bwa mem -A 1 -B 4 -E 50 -L 0 -t $ppn ${index} ${debarcode_path}/${sample}.extract.R2.fq.gz | samtools view -Shb - > ${bam_path}/${sample}_R2.bam

	mkdir $hicMatrix_path
          # $hicFindRestSite --fasta $fa_path --searchPattern GATC -o /date/guomaoni/Reference/HiC_files/GRCm38_MboI_rest_site_positions.bed
          # $hicFindRestSite --fasta $fa_path --searchPattern CT.AG -o /date/guomaoni/Reference/HiC_files/GRCm38_DdeI_rest_site_positions.bed
          # $hicup_digester --re1 ^GATC,MboI:C^TNAG,DdeI --genome GRCm38 --zip --outdir /date/guomaoni/Reference/HiC_files/ $fa_path

	$hicBuildMatrix --samFiles ${bam_path}/${sample}_R1.bam ${bam_path}/${sample}_R2.bam \
	--binSize 10000 \
	--restrictionSequence GATC \
	--danglingSequence GATC \
	--restrictionCutFile ${cutsite} \
	--outBam ${bam_path}/${sample}_hicexplorer.bam \
	--outFileName ${hicMatrix_path}/${sample}_10kb.h5 \
	--QCfolder ${hicMatrix_path}/${sample}_10kb_QC \
	--threads $ppn \
	--inputBufferSize 400000
	 
	### 6.bedtools ###
	mkdir $bedpe_path
	$bedtools bamtobed -bedpe -i ${bam_path}/${sample}_hicexplorer.bam > $bedpe_path/${sample}_hicexplorer.bedpe

	$python_exec $script_path/bedpe2Higashi.py \
	--input_stat=$debarcode_path/raw_stat.csv \
	--input_bed=$bedpe_path/${sample}_hicexplorer.bedpe \
	--output=$bedpe_path/data_all_pixel.txt \
	--stat=$bedpe_path/stat_clean.csv
			 
    ### 7. stat_clean ###
    $r_exec $script_path/stat_clean.R \
	-f $bedpe_path/stat_clean.csv \
	-p $work_path 

        
else
	echo "Please give the mother_path(of the project name) and the project name in order!
	Just like:bash xxx.sh [/media/xxx] [Spatial-E13]
fi

               
