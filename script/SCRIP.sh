# Descargar genoma completo de arabidopsis thaliana

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-61/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O ath_genome.fa.gz
gunzip ath_genome.fa.gz
mv ath_genome.fa.gz referencia.fasta

# Cargar programas necesarios

module load samtools/
module load bwa/
module load bcftools/
 
# Indexar referencia de genoma

samtools faidx referencia.fasta
bwa index referencia.fasta

# Descargar muestra SRA

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR158/085/SRR15810685/SRR15810685_1.fastq.gz -O ath.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1957958/SRR1957958_1.fastq.gz -O cap.fastq.gz

# Realizar filtracion con Trimmimatic

java -jar trimmomatic-0.39.jar SE -phred33 ath.fastq.gz limpio_ath.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
java -jar trimmomatic-0.39.jar SE -phred33 cap.fastq.gz limpio_cap.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Realizar el mapeo con bwa

bwa mem referencia.fasta limpio_ath.fastq.gz > ath.sam
bwa mem referencia.fasta limpio_cap.fastq.gz > cap.sam

# Transformar a Bam y sortear

samtools view -Sb ath.sam | samtools sort -o ath.sorted.bam
samtools view -Sb cap.sam | samtools sort -o cap.sorted.bam

# Indexar

samtools index ath.sorted.bam
samtools index cap.sorted.bam

# Buscar variantes

for sp in ath cap; do
bcftools mpileup -f referencia.fasta ${sp}.sorted.bam | bcftools call -mv -Oz -o ${sp}.vcf.gz
bcftools index ${sp}.vcf.gz
done

# Extraer region genomica del gen FLC

for sp in ath bra cap; do
bcftools mpileup -f referencia.fasta ${sp}.sorted.bam | bcftools call -mv -Oz -o ${sp}.vcf.gz
bcftools index ${sp}.vcf.gz
done

# Generar secuencias consenso en cada especie 

for sp in ath cap; do
bcftools consensus -f FLC_ref.fasta ${sp}.vcf.gz > FLC_${sp}.fasta
done

cat FLC_ath.fasta FLC_bra.fasta FLC_cap.fasta > FLC_todas.fasta

# Alinear con muscle

./muscle3.8.31_i86linux64 -in FLC_todas.fasta -out FLC_alineado.fasta

# Filogenia con iqtree

module load iqtree/2.2.2.6
iqtree2 -s FLC_alineado.fasta

