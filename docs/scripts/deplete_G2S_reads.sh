set -e

BREAKPOINT_BED=$1
GENOME_FASTA=$2
BAM=$3
TELO_SEQ=$4
NUM_THREADS=$5
OUTPUT_DIR=$6


usage() {
    echo "Usage: $0 breakpoint_bed genome_fasta alignment_bam telo_sequence output_dir"
    exit 0
}

if [[ $# != 6 ]]; then usage; fi

TELO_IN_READS=$(printf "${TELO_SEQ}"'%.0s' {1..10})

mkdir -p "${OUTPUT_DIR}"

tmp_genome="${OUTPUT_DIR}/genome.fa"
G2S_bed="${OUTPUT_DIR}/G2S_breakpoints_with_slop.bed"
G2S_reads="${OUTPUT_DIR}/reads_in_G2S_breakpoints.fasta.gz"
telo_reads="${OUTPUT_DIR}/reads__with_soma_specific_telomeres"
non_telo_reads="${OUTPUT_DIR}/reads_without_soma_specific_telomeres"

gzip -dc "${GENOME_FASTA}" > "${tmp_genome}"
samtools faidx "${tmp_genome}"
grep "G2S" "${BREAKPOINT_BED}" | bedtools slop -b 20 -g "${tmp_genome}.fai" > "${G2S_bed}"
# samtools index -@ ${NUM_THREADS} ${BAM}
samtools view -h -@ ${NUM_THREADS} --regions-file "${G2S_bed}" ${BAM} | samtools fasta -0 ${G2S_reads}
seqkit locate -j ${NUM_THREADS} -m 1 -p "${TELO_IN_READS}" ${G2S_reads} | cut -f 1 | tail -n+2 | uniq > ${telo_reads}.txt

# # Telo reads
# samtools view -h -@ ${task.cpus} -N ${telo_reads}.txt ${bam} | samtools fasta -o ${telo_reads}.fasta.gz

## Get no-telo reads (-v flag)

## In BAM format
awk '{print $1"\t1\t"$2}' "${tmp_genome}.fai" > "${tmp_genome}.bed"
samtools view -h -@ ${NUM_THREADS} --regions-file "${tmp_genome}.bed" ${BAM} | samtools view -h -@ ${NUM_THREADS} -N ^${telo_reads}.txt -O BAM -o "${non_telo_reads}.bam"
samtools index -@ ${NUM_THREADS} "${non_telo_reads}.bam"

## In Fasta format
# samtools view -h -@ ${NUM_THREADS} -N ^${telo_reads}.txt ${BAM} | samtools fasta -0 ${non_telo_reads}.fasta.gz

# Remove tmp files
rm "${G2S_reads}" ${tmp_genome}*
