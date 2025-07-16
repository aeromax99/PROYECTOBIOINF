
# Descargar genoma completo de Brassica oleracea

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/416/815/GCA_900416815.2_Brassica_oleracea_HDEM/GCA_900416815.2_Brassica_oleracea_HDEM_genomic.fna.gz -O bra_genome.fna.gz
gunzip bra_genome.fna.gz
mv bra_genome.fna referencia.fasta

# Cargar programas necesarios

source /u/local/Modules/default/init/modules.sh
module load samtools/
module load bwa/
module load bcftools/

trimmo=/u/scratch/d/dechavez/Bioinformatica-PUCE/HerrBio/PROGRAMS/Trimmomatic-0.39/trimmomatic-0.39.jar
cp $trimmo ./
 
# Indexar referencia de genoma

samtools faidx referencia.fasta
bwa index referencia.fasta

# Descargar muestra SRA

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/009/SRR2297049/SRR2297049_1.fastq.gz -O brin.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/001/SRR2297101/SRR2297101_1.fastq.gz -O brol.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/007/SRR2296977/SRR2296977_1.fastq.gz -O brma.fastq.gz

# Realizar filtracion con Trimmimatic

java -jar trimmomatic-0.39.jar SE -phred33 brin.fastq.gz limpio_brin.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
java -jar trimmomatic-0.39.jar SE -phred33 brol.fastq.gz limpio_brol.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
java -jar trimmomatic-0.39.jar SE -phred33 brma.fastq.gz limpio_brma.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Realizar el mapeo con bwa

bwa mem referencia.fasta limpio_brin.fastq.gz > brin.sam
bwa mem referencia.fasta limpio_brol.fastq.gz > brol.sam
bwa mem referencia.fasta limpio_brma.fastq.gz > brma.sam

# Transformar a Bam y sortear

samtools view -Sb brin.sam | samtools sort -o brin.sorted.bam
samtools view -Sb brol.sam | samtools sort -o brol.sorted.bam
samtools view -Sb brma.sam | samtools sort -o brma.sorted.bam

# Indexar

samtools index brin.sorted.bam
samtools index brol.sorted.bam
samtools index brma.sorted.bam

# Buscar variantes

for sp in brin brol brma; do
bcftools mpileup -f referencia.fasta ${sp}.sorted.bam | bcftools call -mv -Oz -o ${sp}.vcf.gz
bcftools index ${sp}.vcf.gz
done

# Descomprimir archivos

gunzip *.vcf.gz

