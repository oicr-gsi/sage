version 1.0

struct GenomeResources {
    String sageModules
    String refFasta
    String ensemblDir
    String hotspots
    String driverGenePanel
    String highConfBed
}

workflow sage {
  input {
    File tumour_bam
    File tumour_bai
    File normal_bam
    File normal_bai
    String donor
    String genomeVersion = "38"
    Boolean use_redux = false
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Int min_map_quality = 10 
    Int hard_min_tumor_qual = 50
    Int hard_min_tumor_raw_alt_support = 2 
    Float hard_min_tumor_vaf = 0.002
  }

  parameter_meta {
    tumour_bam: "Input tumor BAM file"
    tumour_bai: "Input tumor BAI index"
    normal_bam: "Input normal BAM file"
    normal_bai: "Input normal BAI index"
    donor: "Patient/donor identifier"
    genomeVersion: "Genome version (only 38 supported)"
    use_redux: "Run Redux to generate UMI jitter files"
    chromosomes: "List of chromosomes to process in parallel"
  }

  Map[String,GenomeResources] resources = {
    "38": {
      "sageModules": "sage/3.4.4 sage-data/1.0 hg38/p12 hmftools-data/53138 samtools",
      "refFasta": "$HG38_ROOT/hg38_random.fa",
      "ensemblDir": "$HMFTOOLS_DATA_ROOT/ensembl_data",
      "hotspots": "$SAGE_DATA_ROOT/hotspots",
      "driverGenePanel": "$SAGE_DATA_ROOT/driverGenePanel",
      "highConfBed": "$SAGE_DATA_ROOT/highConfidenceBed"
    }
  }

  call extractName as extractTumorName {
    input:
      inputBam = tumour_bam,
      inputBai = tumour_bai,
      modules = "gatk/4.1.6.0 hg38/p12",
      refFasta = resources[genomeVersion].refFasta
  }

  call extractName as extractNormalName {
    input:
      inputBam = normal_bam,
      inputBai = normal_bai,
      modules = "gatk/4.1.6.0 hg38/p12",
      refFasta = resources[genomeVersion].refFasta
  }

  # Run Redux if enabled
  if (use_redux) {
    call redux as reduxTumor {
      input:
        sample_name = extractTumorName.sample_name,
        input_bam = tumour_bam,
        input_bai = tumour_bai,
        refFasta = resources[genomeVersion].refFasta,
        modules = resources[genomeVersion].sageModules
    }

    call redux as reduxNormal {
      input:
        sample_name = extractNormalName.sample_name,
        input_bam = normal_bam,
        input_bai = normal_bai,
        refFasta = resources[genomeVersion].refFasta,
        modules = resources[genomeVersion].sageModules
    }
  }

  # Run SAGE in parallel per chromosome
  scatter (chr in chromosomes) {
    call sagePerChromosome {
      input:
        chromosome = chr,
        tumour_name = extractTumorName.sample_name,
        tumour_bam = tumour_bam,
        tumour_bai = tumour_bai,
        reference_name = extractNormalName.sample_name,
        reference_bam = normal_bam,
        reference_bai = normal_bai,
        refFasta = resources[genomeVersion].refFasta,
        ensemblDir = resources[genomeVersion].ensemblDir,
        hotspots = resources[genomeVersion].hotspots,
        driverGenePanel = resources[genomeVersion].driverGenePanel,
        highConfBed = resources[genomeVersion].highConfBed,
        tumor_jitter = reduxTumor.jitter_file,
        normal_jitter = reduxNormal.jitter_file,
        min_map_quality = min_map_quality,
        hard_min_tumor_qual = hard_min_tumor_qual,
        hard_min_tumor_raw_alt_support = hard_min_tumor_raw_alt_support,
        hard_min_tumor_vaf = hard_min_tumor_vaf,
        modules = resources[genomeVersion].sageModules
    }
  }

  # Merge chromosome VCFs
  call mergeVcfs {
    input:
      vcfs = sagePerChromosome.output_vcf,
      vcf_indices = sagePerChromosome.output_vcf_index,
      sample_name = extractTumorName.sample_name,
      modules = "bcftools/1.9"
  }

  # Merge BQR directories
  call mergeBqrDirs {
    input:
      bqr_zips = sagePerChromosome.output_bqr_directory,
      sample_name = extractTumorName.sample_name
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "SAGE somatic variant calling workflow with Redux UMI processing and per-chromosome parallelization"
    dependencies: [
      {
        name: "SAGE",
        url: "https://github.com/hartwigmedical/hmftools/tree/master/sage"
      },
      {
        name: "Redux",
        url: "https://github.com/hartwigmedical/hmftools/tree/master/redux"
      }
    ]
  }

  output {
    File sage_vcf = mergeVcfs.merged_vcf
    File sage_vcf_index = mergeVcfs.merged_vcf_index
    File sage_bqr_directory = mergeBqrDirs.merged_bqr_zip
    String tumor_sample_name = extractTumorName.sample_name
    String normal_sample_name = extractNormalName.sample_name
    File? tumor_jitter = reduxTumor.jitter_file
    File? normal_jitter = reduxNormal.jitter_file
  }
}

task extractName {
  input {
    File inputBam
    File inputBai
    String modules
    String refFasta
    Int memory = 4
    Int timeout = 4
  }

  command <<<
    set -euo pipefail
    gatk --java-options "-Xmx1g" GetSampleName \
      -R ~{refFasta} \
      -I ~{inputBam} \
      -O sample_name.txt \
      -encode
  >>>

  runtime {
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String sample_name = read_string("sample_name.txt")
  }
}

task redux {
  input {
    String sample_name
    File input_bam
    File input_bai
    String refFasta
    String modules
    Int threads = 8
    Int memory = 16
    Int timeout = 24
  }

  parameter_meta {
    sample_name: "Sample identifier"
    input_bam: "Input BAM file"
    input_bai: "Input BAI index"
    refFasta: "Reference genome FASTA"
    modules: "Required environment modules"
    threads: "Number of threads"
    memory: "Memory in GB"
    timeout: "Timeout in hours"
  }

  command <<<
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
  >>>

  runtime {
    cpu: "~{threads}"
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File jitter_file = "~{sample_name}.redux.jitter.tsv"
  }

  meta {
    output_meta: {
      jitter_file: "Redux UMI jitter statistics file"
    }
  }
}

task sagePerChromosome {
  input {
    String chromosome
    String tumour_name
    File tumour_bam
    File tumour_bai
    String reference_name
    File reference_bam
    File reference_bai
    String refFasta
    String ensemblDir
    String hotspots
    String driverGenePanel
    String highConfBed
    File? tumor_jitter
    File? normal_jitter
    Int min_map_quality
    Int hard_min_tumor_qual
    Int hard_min_tumor_raw_alt_support
    Float hard_min_tumor_vaf
    String modules
    Int threads = 8
    Int memory = 40
    Int timeout = 24
  }

  command <<<
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
  >>>

  runtime {
    cpu: "~{threads}"
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_vcf = "~{tumour_name}.~{chromosome}.sage.vcf.gz"
    File output_vcf_index = "~{tumour_name}.~{chromosome}.sage.vcf.gz.tbi"
    File output_bqr_directory = "~{tumour_name}.~{chromosome}.sage.bqr.zip"
  }
}

task mergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_indices
    String sample_name
    String modules
    Int memory = 8
    Int timeout = 4
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File merged_vcf = "~{sample_name}.sage.vcf.gz"
    File merged_vcf_index = "~{sample_name}.sage.vcf.gz.tbi"
  }
}

task mergeBqrDirs {
  input {
    Array[File] bqr_zips
    String sample_name
    Int memory = 4
    Int timeout = 2
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{memory} GB"
    timeout: "~{timeout}"
  }

  output {
    File merged_bqr_zip = "~{sample_name}.sage.bqr.zip"
  }
}