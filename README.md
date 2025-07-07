# Detección de variantes en el gen FLC en plantas como _Arabidopsis thaliana_ y _Capsella rubella_

### Autor: Andrés Murillo
Email: afmurilloy@puce.edu.ec

## Introduccion 

 El gen FLC es un represor clave de la floración en Arabidopsis thaliana, cuya expresión se modula principalmente por vernalización, un proceso que permite que las plantas florezcan tras la exposición prolongada al frío (Michaels & Amasino, 1999). La regulación del FLC es crucial para sincronizar la floración con condiciones ambientales favorables. Comparar variantes de este gen en especies relacionadas, como Capsella rubella, proporciona información valiosa sobre la evolución de los mecanismos de control del desarrollo floral en plantas Brassicaceae (Slotte et al., 2013). Este tipo de análisis permite identificar polimorfismos funcionales que podrían estar implicados en adaptaciones ecológicas y divergencia fenotípica, siendo especialmente relevante en estudios de genética evolutiva y mejora vegetal. 
 ![FLOR LINDA](https://encrypted-tbn1.gstatic.com/images?q=tbn:ANd9GcSSxRyMsvUdx6_jtw7_N2oZENTbfKp6O24BgdAoReTY6S25F5TIrquINmqTpCD7taycjktcwBUfUS9CdSGWuAxwdg)
 

## Programas requeridos

- SAMTOOLS
-- [`Instalar aquí`](http://www.htslib.org/)
- BWA
-- [`Instalar aquí`](http://bio-bwa.sourceforge.net/)
- BCFTOOLS
-- [`Instalar aquí`](http://www.htslib.org/)
- TRIMMOMATIC
-- [`Instalar aquí`](http://www.usadellab.org/cms/?page=trimmomatic)
- MUSCLE
-- [`Instalar aquí`](https://www.drive5.com/muscle/)
- IQTREE
-- [`Instalar aquí`](http://www.iqtree.org/)

## Program WorkFlow

1. **Descarga e indexación del genoma de referencia**
   - Se descarga el genoma de *Arabidopsis thaliana* (`referencia.fasta`)
   - Indexación con `samtools faidx` y `bwa index`

2. **Descarga de datos de secuenciación**
   - FASTQ de *Arabidopsis thaliana* y *Capsella rubella* desde el SRA

3. **Filtrado de calidad**
   - Uso de `Trimmomatic` para limpiar lecturas: remoción de bases de baja calidad

4. **Mapeo de lecturas**
   - Alineación al genoma con `bwa mem` para obtener archivos `.sam`

5. **Procesamiento de alineamientos**
   - Conversión a `.bam`, ordenamiento e indexación con `samtools`

6. **Llamado de variantes**
   - Con `bcftools mpileup` y `bcftools call` → archivos `.vcf.gz` por especie

7. **Generación de secuencias consenso del gen FLC**
   - A partir de VCFs y una referencia del gen (`FLC_ref.fasta`) usando `bcftools consensus`

8. **Alineamiento múltiple**
   - Se alinean las secuencias consenso con `MUSCLE`

9. **Reconstrucción filogenética**
   - Árbol filogenético generado con `IQ-TREE` usando el alineamiento final


## Instrucciones 

CORRER EL PROGRAMA SCRIPS 

## Resultado esperado

Luego de correr el programa se espera obtener los siguientes archivos: 

- `*.sorted.bam`: lecturas alineadas que serán visibles en IGV 

- `FLC_alineado.fasta.treefile`: árbol filogenético de las especies seleccionadas 


# BIBLIOGRAFÍA 
Michaels, S. D., & Amasino, R. M. (1999). FLOWERING LOCUS C encodes a novel MADS domain protein that acts as a repressor of flowering. The Plant Cell, 11(5), 949–956. https://doi.org/10.1105/tpc.11.5.949

Slotte, T., Hazzouri, K. M., Ågren, J. A., Koenig, D., Maumus, F., Guo, Y. L., ... & Wright, S. I. (2013). The Capsella rubella genome and the genomic consequences of rapid mating system evolution. Nature Genetics, 45(7), 831–835. https://doi.org/10.1038/ng.2669