# sage

SAGE somatic variant calling workflow with Redux UMI processing and per-chromosome parallelization

## Overview

## Dependencies

* [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage)
* [Redux](https://github.com/hartwigmedical/hmftools/tree/master/redux)
* [GATK](https://github.com/broadinstitute/gatk)
* [bcftools](https://github.com/samtools/bcftools)
* [Samtools](https://github.com/samtools/samtools)


## Usage

### Cromwell
```
java -jar cromwell.jar run sage.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumour_bam`|File|Input tumor BAM file
`tumour_bai`|File|Input tumor BAI index
`normal_bam`|File|Input normal BAM file
`normal_bai`|File|Input normal BAI index
`donor`|String|Patient/donor identifier


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`genomeVersion`|String|"hg38"|Genome version (only hg38 supported)
`use_redux`|Boolean|false|Run Redux to generate UMI jitter files
`chromosomes`|Array[String]|["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]|List of chromosomes to process in parallel
`min_map_quality`|Int|10|Minimum map quality
`hard_min_tumor_qual`|Int|50|Minimum hard threshold for tumor base quality
`hard_min_tumor_raw_alt_support`|Int|2|Minimum raw alternate allele support in tumor
`hard_min_tumor_vaf`|Float|0.002|Minimum tumor variant allele frequency


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractTumorName.memory`|Int|4|Memory in GB
`extractTumorName.timeout`|Int|4|Timeout in hours
`extractNormalName.memory`|Int|4|Memory in GB
`extractNormalName.timeout`|Int|4|Timeout in hours
`reduxTumor.threads`|Int|8|Number of threads
`reduxTumor.memory`|Int|16|Memory in GB
`reduxTumor.timeout`|Int|24|Timeout in hours
`reduxNormal.threads`|Int|8|Number of threads
`reduxNormal.memory`|Int|16|Memory in GB
`reduxNormal.timeout`|Int|24|Timeout in hours
`sagePerChromosome.threads`|Int|8|Number of threads
`sagePerChromosome.memory`|Int|40|Memory in GB
`sagePerChromosome.timeout`|Int|24|Timeout in hours
`mergeVcfs.memory`|Int|8|Memory in GB
`mergeVcfs.timeout`|Int|4|Timeout in hours
`mergeBqrDirs.memory`|Int|4|Memory in GB
`mergeBqrDirs.timeout`|Int|2|Timeout in hours


### Outputs

Output | Type | Description | Labels
---|---|---|---
`sage_vcf`|File|Merged VCF file containing somatic variants from all chromosomes|
`sage_vcf_index`|File|Index file for the merged VCF|
`sage_bqr_directory`|File|Merged base quality recalibration directory|
`tumor_jitter`|File?|Optional Redux jitter file for tumor sample|
`normal_jitter`|File?|Optional Redux jitter file for normal sample|


## Commands
 This section lists command(s) run by SAGE workflow
 
 * Running SAGE
 
 
 ```
     set -euo pipefail
     gatk --java-options "-Xmx1g" GetSampleName \
       -R ~{refFasta} \
       -I ~{inputBam} \
       -O sample_name.txt \
       -encode
 ```
 ```
     set -euo pipefail
     
     # Run Redux to generate jitter file
     java -Xmx~{memory}G -jar /.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/workflow/sage/test/redux_v1.2.2.jar \
     -sample ~{sample_name} \
     -input_bam ~{input_bam} \
     -ref_genome ~{refFasta} \
     -ref_genome_version 38 \
     -bamtool $SAMTOOLS_ROOT/bin/samtools \
     -threads ~{threads} \
     -log_level DEBUG \
     -output_dir ./
         
     java -jar redux.jar 
     -sample SAMPLE_ID 
     -input_bam SAMPLE_ID.lane_01.bam,SAMPLE_ID.lane_02.bam,SAMPLE_ID.lane_03.bam  
     -ref_genome /path_to_fasta_files/
     -ref_genome_version V37
     -unmap_regions /ref_data/unmap_regions.37.tsv
     -ref_genome_msi_file /ref_data/msi_jitter_sites.37.tsv.gz 
     -write_stats 
     -bamtool /path_to_samtools/ 
     -output_dir /path_to_output/
     -log_level DEBUG 
     -threads 24
     # Redux outputs: sample.redux.jitter.tsv
 ```
 ```
     set -euo pipefail
     
     mkdir -p ~{tumour_name}.sage.bqr
     
     java -Xmx32G -cp $SAGE_ROOT/sage.jar com.hartwig.hmftools.sage.SageApplication \
       -tumor ~{tumour_name} \
       -tumor_bam ~{tumour_bam} \
       -reference ~{reference_name} \
       -reference_bam ~{reference_bam} \
       -ref_genome_version 38 \
       -ref_genome ~{refFasta} \
       -ensembl_data_dir ~{ensemblDir} \
       -high_confidence_bed ~{highConfBed} \
       ~{if defined(tumor_jitter) then "-tumor_jitter " + tumor_jitter else ""} \
       ~{if defined(normal_jitter) then "-reference_jitter " + normal_jitter else ""} \
       -specific_chr ~{chromosome} \
       -output_vcf ~{tumour_name}.~{chromosome}.sage.vcf.gz \
       -threads ~{threads} \
       -min_map_quality ~{min_map_quality} \
       -hard_min_tumor_qual ~{hard_min_tumor_qual} \
       -hard_min_tumor_raw_alt_support ~{hard_min_tumor_raw_alt_support} \
       -hard_min_tumor_vaf ~{hard_min_tumor_vaf}
 
     # Move BQR files
     mv *.sage.bqr.tsv ~{tumour_name}.sage.bqr/ 2>/dev/null || true
     zip -r ~{tumour_name}.~{chromosome}.sage.bqr.zip ~{tumour_name}.sage.bqr/
 ```
 ```
     set -euo pipefail
     
     # Create file list for bcftools concat
     for vcf in ~{sep=' ' vcfs}; do
       echo "$vcf" >> vcf_list.txt
     done
     
     # Sort by chromosome order
     sort -V vcf_list.txt > vcf_list_sorted.txt
     
     # Concatenate VCFs in order
     bcftools concat \
       --file-list vcf_list_sorted.txt \
       --output-type z \
       --output ~{sample_name}.sage.vcf.gz
     
     # Index the merged VCF
     tabix -p vcf ~{sample_name}.sage.vcf.gz
 ```
 ```
     set -euo pipefail
     
     mkdir -p ~{sample_name}.sage.bqr
     
     # Unzip all BQR files into the same directory
     for bqr_zip in ~{sep=' ' bqr_zips}; do
       unzip -o "$bqr_zip" -d temp_bqr/
     done
     
     # Move all .tsv files to final directory
     find temp_bqr -name "*.sage.bqr.tsv" -exec mv {} ~{sample_name}.sage.bqr/ \;
     
     # Zip the merged directory
     zip -r ~{sample_name}.sage.bqr.zip ~{sample_name}.sage.bqr/
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
