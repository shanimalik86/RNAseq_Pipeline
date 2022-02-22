configfile: "config.yaml"

envvars:
	"SLURM_JOB_ID"

DATA_DIR = config["data"]

if config["trim"]["type"]=="SE":
	SAMPLES, = glob_wildcards(DATA_DIR + "/{sample}.fastq.gz")
elif config["trim"]["type"]=="PE":
	SAMPLES, = glob_wildcards(DATA_DIR + "/{sample}_R1.fastq.gz")
else:
	raise ValueError('please specify only "SE" or "PE" for the "strand" parameter in the config file.')
	
rule all:
	input:
		counting = expand("counts/{sample}_featcounts.txt", sample=SAMPLES)
	message:
		"i am all here"
        
if config["trim"]["type"]=="SE":
	rule trim:
		input: 
			DATA_DIR + "/{sample}.fastq.gz"
		output:
			"trimmomatic_fastq/{sample}.trimmed.fastq.gz"
		message:
			"I am trimming"
		params:
			length=config["trim"]["length"],
			trail=config["trim"]["trail"],
			lead=config["trim"]["lead"],
			window=config["trim"]["window"],
			strand=config["trim"]["strand"]
		shell:
			"""
			mkdir -p trimmomatic_fastq
			module load Trimmomatic
			trimmomatic {params.strand} -threads 2 -phred33 {input} {output} ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:{params.window} LEADING:{params.lead} TRAILING:{params.trail} MINLEN:{params.length}
			"""
else:
	rule trim:
		input:
			r1=DATA_DIR + "/{sample}_R1.fastq.gz",
			r2=DATA_DIR + "/{sample}_R2.fastq.gz"
		output:
			TR1="trimmomatic_fastq/{sample}_R1.trimmed.fastq.gz",
			TR2="trimmomatic_fastq/{sample}_R2.trimmed.fastq.gz",
			TR1un="trimmomatic_fastq/{sample}_R1un.trimmed.fastq.gz",
			TR2un="trimmomatic_fastq/{sample}_R2un.trimmed.fastq.gz"
		message:
			"I am paired trimming"
		params:
			length=config["trim"]["length"],
                        trail=config["trim"]["trail"],
                        lead=config["trim"]["lead"],
                        window=config["trim"]["window"],
                        strand=config["trim"]["type"]
		shell:
			"""
			mkdir -p trimmomatic_fastq
                        module load Trimmomatic
                        trimmomatic {params.strand} -threads 2 -phred33 {input.r1} {input.r2} {output.TR1} {output.TR1un} {output.TR2} {output.TR2un} ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:{params.window} LEADING:{params.lead} TRAILING:{params.trail} MINLEN:{params.length}
                        """				
rule STAR_index:	
	input:
		fasta = config["fasta"]
	output:
		genome_index = ["genome/"+ f for f in ["chrLength.txt","chrNameLength.txt","chrName.txt","chrStart.txt","Genome","genomeParameters.txt","SA","SAindex"]]
	
	params:
		genome_dir = "genome/",
		index_id = os.environ["SLURM_JOB_ID"]
	
	message:
		"I am STAR index"
	
	threads: 10

	shell:
		"""
		mkdir -p {params.genome_dir}
		module load STAR/2.7.3a-IGB-gcc-8.2.0
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genome_dir} --genomeFastaFiles {input.fasta} --limitGenomeGenerateRAM 32000000000 --outTmpDir /scratch/{params.index_id}
		"""

if config["trim"]["type"]=="SE":

	rule mapping:
		input:
			trimfile = "trimmomatic_fastq/{sample}.trimmed.fastq.gz",
			genome_index = rules.STAR_index.output
		output:
			"star/{sample}_Aligned.sortedByCoord.out.bam",
			"star/{sample}_Log.final.out"
		params:
			prefix= "star/{sample}_",
			genome_index = "genome/",
			gtf = config["gtf"],
			job_id= os.environ["SLURM_JOB_ID"]
	
		threads: 10
    
		shell:
			"""
			mkdir -p star
			module load STAR/2.7.3a-IGB-gcc-8.2.0
			STAR --runThreadN {threads} --readFilesCommand zcat --genomeDir {params.genome_index} --readFilesIn {input.trimfile} --sjdbGTFfile {params.gtf} --sjdbOverhang 99 --outFileNamePrefix {params.prefix} --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outTmpDir /scratch/{params.job_id}
			"""
else:
	rule mapping:
                input:
                        TR1 = "trimmomatic_fastq/{sample}_R1.trimmed.fastq.gz",
                        TR2 = "trimmomatic_fastq/{sample}_R2.trimmed.fastq.gz",
                        genome_index = rules.STAR_index.output
                output:
                        "star/{sample}_Aligned.sortedByCoord.out.bam",
                        "star/{sample}_Log.final.out"
                params:
                        prefix= "star/{sample}_",
                        genome_index = "genome/",
                        gtf = config["gtf"],
                        job_id= os.environ["SLURM_JOB_ID"]

                message:
                        "i am paired mapping"

                threads: 10

                shell:
                	"""
                	mkdir -p star
                	module load STAR/2.7.3a-IGB-gcc-8.2.0
                	STAR --runThreadN {threads} --readFilesCommand zcat --genomeDir {params.genome_index} --readFilesIn {input.TR1} {input.TR2} --sjdbGTFfile {params.gtf} --sjdbOverhang 99 --outFileNamePrefix {params.prefix} --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outTmpDir /scratch/{params.job_id}
			"""
rule counts:
	input:
		bams = "star/{sample}_Aligned.sortedByCoord.out.bam",
		gtf = config["gtf"]
	output:
		counting="counts/{sample}_featcounts.txt"
	message:
		"I am counting"
	shell:
		"""
		mkdir -p counts
		module load Subread
		featureCounts -T 4 -s 2 -g gene_id -t exon -o {output.counting} -a {input.gtf} {input.bams}
		"""
