#!/bin/bash

######## Step 1: Check files and parse command line arguments #######
## Set default values (Please change to your own diretories):
WORKDIR="/home/scratch/dn_bulk"
REFDIR="/home/users/haojiawu/refdata-cellranger-GRCh38-3.0.0"
FASTQDIR="/home/scratch/dn_bulk/fastq_dir"
do_fastqc=false
paired=false
n_cores=20
## Set command line arguments
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		-g|--genome)
		REFDIR="$2"
		shift # past argument
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift # past argument
		;;
                -f|--fqdir)
                FASTQDIR="$2"
                shift # past argument
                ;;
		-t|--cpus)
		n_cores="$2"
		shift # past argument
		;;        
		--do_fastqc)
		do_fastqc=true
		shift # past argument
		;;
		-p|--paired)
		paired=false
		shift # past argument
		;;
		-h|--help)
		echo "Usage: ./run_RNAseq.sh -g <REFDIR> -w <WORKDIR> -f <FASTQDIR>"
		echo "Optional:"
		echo "  --do_fastqc: FastQC to check the read quality. Default is false"
		echo "  -t, --cpus: Number of cores to use. Default is 12"
		exit
		;;
		*)
		# unknown option
		echo "Unknown option: $key, exiting."
		echo "Usage: ./run_RNAseq.sh -g <REFDIR> -w <WORKDIR> -f <FASTQDIR>"
		echo "Optional:"
		echo "  --do_fastqc: FastQC to check the read quality. Default is false"
		echo "  -t, --cpus: Number of cores to use. Default is 12"
		exit
		;;
	esac
	shift
done

## Check the working directory
if [[ ! -d $WORKDIR ]]; then
echo "The working directory does not exist: $WORKDIR. Exiting..."
	exit 1
elif [[ ! -d $REFDIR ]]; then
echo "The reference directory does not exist: $REFDIR. Exiting..."
        exit 2
elif [[ ! -d $FASTQDIR ]]; then
echo "The fastq directory does not exist: $FASTQDIR. Exiting..."
        exit 3
else
	## if the drectory exists, turn it to a absolute directory
	REFDIR=$(readlink -e $REFDIR)
	WORKDIR=$(readlink -e $WORKDIR)
        FASTQDIR=$(readlink -e $FASTQDIR)
	echo "REFDIR=$REFDIR"
	echo "WORKDIR=$WORKDIR"
    echo "FASTQDIR=$FASTQDIR"
fi

## Check the fasta and gtf files
if [[ ! -d $REFDIR ]]; then
	echo "The reference directory does not exist: $REFDIR. Exiting..."
	exit 1
else
	REF_GTF="$REFDIR/genes/genes.gtf"
	REF_FASTA="$REFDIR/fasta/genome.fa"
	if [[ ! -f $REF_GTF ]]; then
		echo "$REF_GTF file not found. Please include a gtf file in the reference directory."
		exit 1
	fi
	if [[ ! -f $REF_FASTA ]]; then
		echo "$REF_FASTA file not found. Please include a fasta file in the reference directory."
		exit 1
	fi
	STAR_INDEX="$REFDIR/star/"
fi


cd $WORKDIR

## Create the output directories
mkdir -p fastQC_dir
mkdir -p trimgalore_dir
mkdir -p star_dir
mkdir -p featureCounts_dir/exon_dir
mkdir -p featureCounts_dir/intron_dir


########## Step 2: Read QC, trimming, and alignment ##########
echo "Read QC and trimming..."
if [ paired = false ]; then
	cd fastq_dir
	for fq in $(ls); do
		if [ "$do_fastqc" = false ]; then
			echo "Performing FastQC for $fqname"
			fastqc $fq -o ../fastQC_dir
            trim_galore --fastqc --fastqc_args "--outdir ../fastQC_dir" --output_dir ../trimgalore_dir  $fq
        else
        trim_galore --output_dir ../trimgalore_dir $fq
		fi
        fqname=$(echo $fq | cut -f1 -d '.')
        append="_trimmed.fq.gz"
        suffix="/home/scratch/dn_bulk/trimgalore_dir/"
        fq_trimmed="$suffix$fqname$append"

echo "Mapping reads from $fqname to the reference genome mm10..."
	STAR \
	    --genomeDir $STAR_INDEX \
	    --sjdbGTFfile $REF_GTF \
	    --runThreadN $n_cores \
	    --outFileNamePrefix ../star_dir/$fqname. \
	    --readFilesIn $fq_trimmed \
            --sjdbOverhang 100 \
            --readFilesCommand zcat \
            --outSAMmultNmax 1 \
            --outFilterMultimapNmax 50 \
            --outSAMunmapped Within \
            --twopassMode Basic 
   cd ../star_dir
   samtools view -@ $n_cores -S -b $fqname.Aligned.out.sam -o $fqname.filtered.tagged.Aligned.out.bam
   rm $fqname.Aligned.out.sam
   cd ../
   Rscript --vanilla Featurecountparameters.R -g $REF_GTF -b $WORKDIR/star_dir/$fqname.filtered.tagged.Aligned.out.bam -c $n_cores
   cd star_dir
   mv $fqname.ex.stat.txt ../featureCounts_dir/exon_dir/
   mv $fqname.in.stat.txt ../featureCounts_dir/intron_dir/
   mv $fqname.ex.counts.txt ../featureCounts_dir/exon_dir/
   mv $fqname.in.counts.txt ../featureCounts_dir/intron_dir/
   rm *.bam
   cd ../trimgalore_dir
   rm *.fq.gz
   cd ../fastq_dir
   done

cd $WORKDIR
   Rscript Count_matrix.R

fi

if [ paired = true ]; then

	cd fastq_dir
  for basename in $(ls | cut -f1 -d '_' | sort | uniq); do
		echo $basename
		fq1="_1.fastq.gz"
		fq2="_2.fastq.gz"
		fq1=$basename$fq1
		fq2=$basename$fq2
		if [ "$do_fastqc" = false ]; then
			echo "Performing FastQC for $basename"
			fastqc $fq1 -o ../fastQC_dir
			fastqc $fq2 -o ../fastQC_dir
           trim_galore --fastqc --fastqc_args "--outdir ../fastQC_dir" --paired --output_dir ../trimgalore_dir  $fq1 $fq2
        else
           trim_galore --paired --output_dir ../trimgalore_dir $fq1 $fq2
		fi
        append1="_1_val_1.fq.gz"
        append2="_2_val_2.fq.gz"
        suffix="/home/scratch/dn_bulk/trimgalore_dir/"
        fq_trimmed1="$suffix$basename$append1"
        fq_trimmed2="$suffix$basename$append2"

echo "Mapping reads from $fqname to the reference genome hg19..."
	/home/users/haojiawu/cellranger-3.1.0/STAR/5dda596/STAR \
	    --genomeDir $STAR_INDEX \
	    --sjdbGTFfile $REF_GTF \
	    --runThreadN $n_cores \
	    --outFileNamePrefix ../star_dir/$basename. \
	    --readFilesIn $fq_trimmed1 $fq_trimmed2  \
      --sjdbOverhang 100 \
      --readFilesCommand zcat \
      --outSAMmultNmax 1 \
      --outFilterMultimapNmax 50 \
      --outSAMunmapped Within \
      --twopassMode Basic 
   cd ../star_dir
   samtools view -@ $n_cores -S -b $basename.Aligned.out.sam -o $basename.filtered.tagged.Aligned.out.bam
   rm $basename.Aligned.out.sam
   cd ../
   Rscript --vanilla Featurecountparameters.R -g $REF_GTF -b $WORKDIR/star_dir/$basename.filtered.tagged.Aligned.out.bam -c $n_cores
   cd star_dir
   mv $basename.ex.stat.txt ../featureCounts_dir/exon_dir/
   mv $basename.in.stat.txt ../featureCounts_dir/intron_dir/
   mv $basename.ex.counts.txt ../featureCounts_dir/exon_dir/
   mv $basename.in.counts.txt ../featureCounts_dir/intron_dir/
   rm *.bam
   cd ../trimgalore_dir
   rm *.fq.gz
   cd ../fastq_dir
   done

cd $WORKDIR
   Rscript Count_matrix.R

fi
