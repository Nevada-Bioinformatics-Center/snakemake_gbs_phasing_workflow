import glob

configfile: "config.yaml"

inputbam=config["input_data"]
ref=config["ref"]
freebayesparams=config["freebayes"]

##### target rules #####
rule all:
    input: 
       "calls/all_merged.vcf",
       "calls/all_merged.vcf.stats",
       "calls/all_merged.g70mac3.recode.vcf.stats",
       "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.stats",
       "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf.stats",
       "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
       "calls/all_merged.g70mac3lm.recode.vcf.stats",
       "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.gz",
       "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf",
       "plink/fullset/plink.genome",
       "plink/noindels/plink.genome",
       "plink/myplink.ped",
       "plink/noindels/ssr_progeny_population_noindels.chr-1.ped",
       "plink/fullset/ssr_progeny_population.chr-1.ped",
       "plink/fullset/plink.blocks",
       "plink/fullset/plink.ld",
       "plink/fullset/plink.prune.in",
       "plink/noindels/plink.blocks",
       "plink/noindels/plink.ld",
       "plink/noindels/plink.prune.in",

       ## OLD - can uncomment to look at other filtering parameters though
       #"calls/all_merged.g50mac3.recode.vcf",
       #"calls/all_merged.g50mac3lm.recode.vcf",
       #"calls/all_merged.g25mac3.recode.vcf",
       #"calls/all_merged.g75mac3.recode.vcf",
       #"calls/all_merged.g90mac3.recode.vcf",
       #"calls/all_merged.g10mac3.recode.vcf",
       #"calls/all_merged.g100mac3.recode.vcf",
       #"calls/all_merged.g90mac3lm.recode.vcf",
       #"calls/all_merged.g50mac3lm.recode.fixed.phased.vcf.gz",
       #"calls/all_merged.g90mac3lm.recode.fixed.phased.vcf.gz",
       #"calls/all_merged.g50mac3lm.noindel.recode.vcf",
       #"calls/all_merged.g50mac3lm.noindel.recode.fixed.vcf",
       #"calls/all_merged.g90mac3lm.fixed.refremoved.recode.vcf",
       #"calls/all_merged.g50mac3lm.fixed.refremoved.recode.vcf",
       #"calls/all_merged.g50mac3lm.recode.phased.imputed.vcf.gz",
       #"calls/all_merged.unobs.vcf",
       #"calls/all_merged.unobs.vcf.stats",


rule sambamba_markdup:
    input:
        inputbam
    output:
        "mapped/all_merged.rmdup.bam"
    params:
        extra="-r"  # optional parameters
    log:
        "logs/sambamba-markdup/all_merged.log"
    threads: 16
    resources: time_min=480, mem_mb=40000, cpus=16, tmpdir="tmp/"
    wrapper:
        "0.73.0/bio/sambamba/markdup"


rule samtools_index_merged:
    input:
        "mapped/all_merged.rmdup.bam"
    output:
        "mapped/all_merged.rmdup.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.73.0/bio/samtools/index"

rule freebayes:
    input:
        ref=ref,
        samples="mapped/all_merged.rmdup.bam",
        indexes="mapped/all_merged.rmdup.bam.bai"
    output:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/all_merged.log"
    params:
        extra=freebayesparams,         # optional parameters
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # flag to use bcftools norm to normalize indels
    threads: 16
    resources: time_min=480, mem_mb=40000, cpus=16, tmpdir="tmp/"
    wrapper:
        "0.73.0/bio/freebayes"

rule freebayes_unobs:
    input:
        ref=ref,
        samples="mapped/all_merged.rmdup.bam",
        indexes="mapped/all_merged.rmdup.bam.bai"
    output:
        "calls/all_merged.unobs.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/all_merged_unobs.log"
    params:
        extra="--haplotype-length 3 --min-base-quality 10 --read-mismatch-limit 3 --min-coverage 5 --exclude-unobserved-genotypes --genotype-qualities --mismatch-base-quality-threshold 10 -u",         # optional parameters
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # flag to use bcftools norm to normalize indels
    threads: 16
    resources: time_min=480, mem_mb=40000, cpus=16, tmpdir="tmp/"
    wrapper:
        "0.73.0/bio/freebayes"


rule bcftools_stats:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"


rule bcftools_stats_unobs:
    input:
        "calls/all_merged.unobs.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.unobs.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule vcftools_filter_maxmiss:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g50mac3.recode.vcf"
    params:
        "calls/all_merged.g50mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss50.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.5 --mac 3 --minQ 20 --maf 0.03 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss75:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g75mac3.recode.vcf"
    params:
        "calls/all_merged.g75mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss75.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.75 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss70:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3.recode.vcf"
    params:
        "calls/all_merged.g70mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss70.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.70 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss25:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g25mac3.recode.vcf"
    params:
        "calls/all_merged.g25mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss25.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.25 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss90:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g90mac3.recode.vcf"
    params:
        "calls/all_merged.g90mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss90.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.90 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss10:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g10mac3.recode.vcf"
    params:
        "calls/all_merged.g10mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss10.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 0.10 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_maxmiss100:
    input:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g100mac3.recode.vcf"
    params:
        "calls/all_merged.g100mac3"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/maxmiss100.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --max-missing 1 --mac 3 --minQ 20 --maf 0.05 --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_missind:
    input:
        "calls/all_merged.g50mac3.recode.vcf"
    output:
        "calls/all_merged.g50mac3.imiss"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/missindv.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --missing-indv --stdout > {output} 2> {log}"

rule vcftools_missind90:
    input:
        "calls/all_merged.g90mac3.recode.vcf"
    output:
        "calls/all_merged.g90mac3.imiss"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/missindv90.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --missing-indv --stdout > {output} 2> {log}"

rule vcftools_missind70:
    input:
        "calls/all_merged.g70mac3.recode.vcf"
    output:
        "calls/all_merged.g70mac3.imiss"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/missindv70.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --missing-indv --stdout > {output} 2> {log}"

rule mawk_missind:
    input:
        "calls/all_merged.g50mac3.imiss"
    output:
        "calls/lowDP50.indv"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/vcftools.yaml"
    shell:
        "mawk '$5 > 0.5' {input} | cut -f1 > {output}"

rule mawk_missind90:
    input:
        "calls/all_merged.g90mac3.imiss"
    output:
        "calls/lowDP90.indv"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/vcftools.yaml"
    shell:
        "mawk '$5 > 0.5' {input} | cut -f1 > {output}"

rule mawk_missind70:
    input:
        "calls/all_merged.g70mac3.imiss"
    output:
        "calls/lowDP70.indv"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/vcftools.yaml"
    shell:
        "mawk '$5 > 0.5' {input} | cut -f1 > {output}"

rule vcftools_filter_missind:
    input:
        vcf="calls/all_merged.g50mac3.recode.vcf",
	miss="calls/lowDP50.indv"
    output:
        "calls/all_merged.g50mac3lm.recode.vcf"
    params:
        "calls/all_merged.g50mac3lm"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/removemissindv.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.miss} --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_noindels50:
    input:
        "calls/all_merged.g50mac3lm.recode.vcf"
    output:
        "calls/all_merged.g50mac3lm.noindel.recode.vcf"
    params:
        "calls/all_merged.g50mac3lm.noindel"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/remove_indels50.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input} --remove-indels --recode --recode-INFO-all --out {params} 2> {log}"

rule fix_missind50_noindel:
    input:
        "calls/all_merged.g50mac3lm.noindel.recode.vcf"
    output:
        "calls/all_merged.g50mac3lm.noindel.recode.fixed.vcf"
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "grep -v super {input} |  perl -pe 's/\s\.:/\t.\/.:/g'  > {output}"

rule fix_missind50:
    input:
        "calls/all_merged.g50mac3lm.recode.vcf"
    output:
        "calls/all_merged.g50mac3lm.recode.fixed.vcf"
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "grep -v super {input} |  perl -pe 's/\s\.:/\t.\/.:/g'  > {output}"

rule vcftools_filter_missind90:
    input:
        vcf="calls/all_merged.g90mac3.recode.vcf",
	miss="calls/lowDP90.indv"
    output:
        "calls/all_merged.g90mac3lm.recode.vcf"
    params:
        "calls/all_merged.g90mac3lm"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/removemissindv90.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.miss} --recode --recode-INFO-all --out {params} 2> {log}"

rule vcftools_filter_missind70:
    input:
        vcf="calls/all_merged.g70mac3.recode.vcf",
	miss="calls/lowDP70.indv"
    output:
        "calls/all_merged.g70mac3lm.recode.vcf"
    params:
        "calls/all_merged.g70mac3lm"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/removemissindv70.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.miss} --recode --recode-INFO-all --out {params} 2> {log}"


##Beagle has a problem recognizing . as ./. unknown genotype.  I had to recode
## Also beagle has a problem when there aren't more than 1 SNP on a given contig to phase. This is largely a problem
## with the Super contigs such as Super_18. We will remove these
rule fix_missind90:
    input:
        "calls/all_merged.g90mac3lm.recode.vcf"
    output:
        "calls/all_merged.g90mac3lm.recode.fixed.vcf"
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "grep -v super {input} |  perl -pe 's/\s\.:/\t.\/.:/g'  > {output}"

rule beagle_phase_missind90:
    input:
        "calls/all_merged.g90mac3lm.recode.fixed.vcf"
    output:
        "calls/all_merged.g90mac3lm.recode.fixed.phased.vcf.gz"
    params:
        "calls/all_merged.g90mac3lm.recode.fixed.phased"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/beagle/phase_missindv90.log"
    shell:
        "java -jar beagle/beagle.28Jun21.220.jar gt={input} out={params} > {log} 2>&1"

rule fix_missind70:
    input:
        "calls/all_merged.g70mac3lm.recode.vcf"
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.vcf"
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "grep -v super {input} |  perl -pe 's/\s\.:/\t.\/.:/g'  > {output}"

rule beagle_phase_missind70:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.vcf"
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.gz"
    params:
        "calls/all_merged.g70mac3lm.recode.fixed.phased"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/beagle/phase_missindv70.log"
    shell:
        "java -jar beagle/beagle.28Jun21.220.jar gt={input} out={params} > {log} 2>&1"

rule vcftools_filter_noindels70_afterphase:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.gz"
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf"
    params:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/remove_indels70.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --remove-indels --recode --recode-INFO-all --out {params} 2> {log}"


##50
rule beagle_phase_missind50:
    input:
        "calls/all_merged.g50mac3lm.recode.fixed.vcf"
    output:
        "calls/all_merged.g50mac3lm.recode.fixed.phased.vcf.gz"
    params:
        "calls/all_merged.g50mac3lm.recode.fixed.phased"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/beagle/phase_missindv50.log"
    shell:
        "java -jar beagle/beagle.28Jun21.220.jar gt={input} out={params} > {log} 2>&1"


rule vcftools_missind90_afterphase:
    input:
        "calls/all_merged.g90mac3lm.recode.fixed.phased.vcf.gz",
    output:
        "calls/all_merged.g90mac3lm.recode.fixed.phased.imiss"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/missindv90_afterphase.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --stdout > {output} 2> {log}"

rule vcftools_missind50_afterphase:
    input:
        "calls/all_merged.g50mac3lm.recode.fixed.phased.vcf.gz",
    output:
        "calls/all_merged.g50mac3lm.recode.fixed.phased.imiss"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/missindv50_afterphase.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --stdout > {output} 2> {log}"

rule mawk_allind90:
    input:
        "calls/all_merged.g90mac3lm.recode.fixed.phased.imiss"
    output:
        "calls/all90.indv"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/vcftools.yaml"
    shell:
        "cut -f1 {input} > {output}"

rule mawk_allind50:
    input:
        "calls/all_merged.g50mac3lm.recode.fixed.phased.imiss"
    output:
        "calls/all50.indv"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/vcftools.yaml"
    shell:
        "cut -f1 {input} > {output}"


rule vcftools_filter_missind50_afterphase:
    input:
        vcf="calls/all_merged.g50mac3lm.recode.fixed.vcf",
	miss="calls/all50.indv"
    output:
        "calls/all_merged.g50mac3lm.fixed.refremoved.recode.vcf"
    params:
        "calls/all_merged.g50mac3lm.fixed.refremoved"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/remove_ind50.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.miss} --recode --recode-INFO-all --out {params} 2> {log}"


rule vcftools_filter_missind90_afterphase:
    input:
        vcf="calls/all_merged.g90mac3lm.recode.fixed.vcf",
	miss="calls/all90.indv"
    output:
        "calls/all_merged.g90mac3lm.fixed.refremoved.recode.vcf"
    params:
        "calls/all_merged.g90mac3lm.fixed.refremoved"
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/vcftools/remove_ind90.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.miss} --recode --recode-INFO-all --out {params} 2> {log}"



rule bcftools_stats_g70mac3:
    input:
        "calls/all_merged.g70mac3.recode.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3.recode.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_stats_g70mac3lm:
    input:
        "calls/all_merged.g70mac3lm.recode.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3lm.recode.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_stats_g70mac3lmphased:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.gz"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_stats_g70mac3lmphased_noindels:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf.stats",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_g70mac3lmphased_namedvariants:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.vcf.gz"  # either .vcf or .bcf
    output:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
    priority: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/bcftools/namedvariants_70.log"
    conda:
        "envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input} && bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {input} > {output} 2> {log}"



##PLINK section
rule plink_mkdir:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
    output:
        directory("plink/")
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "mkdir -p {input}"

rule plink_mkdir_noindels:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
    output:
        directory("plink/noindels/")
    resources: time_min=320, mem_mb=8000, cpus=1
    shell:
        "mkdir -p {input}"

rule plink_dst70:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
    output:
        "plink/fullset/plink.genome",
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/plink/plink_dst70.log",
    conda:
        "envs/plink.yaml"
    shell:
        "plink --vcf {input} --genome 2> {log} && mv plink.genome plink/fullset/"


rule plink_dst70_noindels:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.noindels.recode.vcf",
    output:
        "plink/noindels/plink.genome",
    resources: time_min=320, mem_mb=8000, cpus=1
    log:
        "logs/plink/plink_dst70_noindels.log",
    conda:
        "envs/plink.yaml"
    shell:
        "plink --vcf {input[0]} --genome 2> {log} && mv plink.genome plink/noindels/"

rule plink_converto_plink70:
    input:
        "calls/all_merged.g70mac3lm.recode.fixed.phased.named.vcf",
    output:
        "plink/myplink.ped",
        "plink/myplink.map",
    resources: time_min=320, mem_mb=8000, cpus=5
    log:
        "logs/plink/plink_convertto_plink70.log",
    conda:
        "envs/plink.yaml"
    shell:
        "plink --threads 5 --vcf {input} --double-id --recode --out myplink 2> {log} && mv myplink* plink/"

rule plink_converto_linkage70:
    input:
        "plink/myplink.ped",
    output:
        "plink/fullset/ssr_progeny_population.chr-1.ped",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --recode HV --out ssr_progeny_population && mv ssr_progeny_population* fullset/"

rule plink_converto_linkage70_noindels:
    input:
        "plink/myplink.ped",
    output:
        "plink/noindels/ssr_progeny_population_noindels.chr-1.ped",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --recode HV --snps-only just-acgt --out ssr_progeny_population_noindels  && mv ssr_progeny_population_noindels* noindels/"

rule plink_indeppairwise70:
    input:
        "plink/myplink.ped",
    output:
        "plink/fullset/plink.prune.in",
        "plink/fullset/plink.prune.out",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --indep-pairphase 50 5 0.5 && mv plink.prune* fullset/"

rule plink_ld70:
    input:
        "plink/myplink.ped",
    output:
        "plink/fullset/plink.ld",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --r2 --ld-window-r2 0 && mv plink.ld fullset/"

rule plink_blocks70:
    input:
        "plink/myplink.ped",
    output:
        "plink/fullset/plink.blocks",
        "plink/fullset/plink.blocks.det",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --blocks no-pheno-req && mv plink.blocks* fullset/"

rule plink_indeppairwise70_noindels:
    input:
        "plink/myplink.ped",
    output:
        "plink/noindels/plink.prune.in",
        "plink/noindels/plink.prune.out",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --indep-pairphase 50 5 0.5 && mv plink.prune* noindels/"

rule plink_ld70_noindels:
    input:
        "plink/myplink.ped",
    output:
        "plink/noindels/plink.ld",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --r2 --ld-window-r2 0 && mv plink.ld noindels/"

rule plink_blocks70_noindels:
    input:
        "plink/myplink.ped",
    output:
        "plink/noindels/plink.blocks",
        "plink/noindels/plink.blocks.det",
    params:
        "myplink"
    resources: time_min=320, mem_mb=8000, cpus=1
    conda:
        "envs/plink.yaml"
    shell:
        "cd plink/ && plink --file {params} --blocks no-pheno-req && mv plink.blocks* noindels/"
