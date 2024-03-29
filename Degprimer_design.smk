import os

configfile: "./config_primer_design.json"

ds=config["degeneracies"].split(",")
ls=config["primer_sizes"].split(",")

rule all:
  input: expand("Results/List_{suffix}.txt", suffix=config["Output_suffix_name"])

rule gff_file:
  input: config["tar_gff_files"]
  output: config["gff_file"]
  params: n=config["NCBI_gene_name"], t=os.path.join(config["work_dir"], "temporal_gff"), o=config["gff_file"], r=config["NCBI_gene_region"]
  message: "extracting positions from annotation file"
  shell: "bash src/create_gff.sh {input} {params.t} {params.n} {params.r} {params.o}"

rule gene_sequences:
  input: gff=config["gff_file"], g=config["tar_file_genomes"]
  output: config["gene_sequences"]
  params: t=os.path.join(config["work_dir"], "temporal_genomes"), o=config["gene_sequences"], n=config["NCBI_gene_name"]
  message: "Extracting gene sequences from genomes"
  shell:  '''
              bash src/extract_gff.sh {params.t} {input.g}
              python src/from_refseq_gff_to_fasta.py -i {input.gff} -d {params.t} -o {params.o} -n {params.n}
              rm -r {params.t}
          '''
rule msa:
  input: config["gene_sequences"]
  output: os.path.join(config["work_dir"], "MSA", "_".join(["MSA_", config["Output_suffix_name"]]))
  params: config["muscle_params"]
  message: "Multiple sequence alignment"
  shell: "muscle -in {input} -out {output} {params}"

rule degeprime_trim:
  input: os.path.join(config["work_dir"], "MSA", "_".join(["MSA_", config["Output_suffix_name"]]))
  output: os.path.join(config["work_dir"], "MSA", "_".join(["trimmed_align_file_from_MSA", config["Output_suffix_name"]]))
  params: config["trim_degeprime"]
  shell: "perl src/DEGEPRIME/TrimAlignment.pl -i {input} {params} -o {output}"

rule multidegeprime:
  input: os.path.join(config["work_dir"], "MSA", "_".join(["trimmed_align_file_from_MSA", config["Output_suffix_name"]]))
  output: all=expand("degeprime_output/output_file_from_MSA_{suffix}_uniq_d{{d}}_l{{l}}", suffix=config["Output_suffix_name"], d=ds, l=ls),
          sl=expand("degeprime_output/Selected_{suffix}/output_file_from_MSA_{suffix}_uniq_d{{d}}_l{{l}}_selected", suffix=config["Output_suffix_name"], d=ds, l=ls)
  params: of=os.path.join(config["work_dir"], "degeprime_output", "output_file_from_MSA_"+config["Output_suffix_name"]+"_uniq_d{d}_l{l}"),
          osf="output_file_from_MSA_"+config["Output_suffix_name"]+"_uniq_d{d}_l{l}_selected",
          osd=os.path.join(config["work_dir"], "degeprime_output", "Selected_"+config["Output_suffix_name"]),
          log=os.path.join(config["work_dir"], "degeprime_output", "Log_files","output_file_from_MSA_"+config["Output_suffix_name"]+"_uniq_d{d}_l{l}.log"),
          ds="{d}", ls="{l}"
  log: expand("degeprime_output/Log_files/output_file_from_MSA_{suffix}_uniq_d{{d}}_l{{l}}.log", suffix=config["Output_suffix_name"], d=ds, l=ls)
  message: "DEGEPRIME -Finding primers per position -d {params.ds} -l {params.ls} "
  shell: '''
            perl src/DEGEPRIME/DegePrime.pl -i {input} -d {params.ds} -l {params.ls} -o {params.of} > {params.log}
            python src/parse_degeprime.py -i {params.of} -o {params.osf} -d {params.osd}
         '''

rule primer_screen_list:
  input: expand("degeprime_output/Selected_{suffix}/output_file_from_MSA_{suffix}_uniq_d{d}_l{l}_selected", suffix=config["Output_suffix_name"], d=ds, l=ls)
  output: expand("potential_primers/primers_{suffix}", suffix=config["Output_suffix_name"])
  params: dir=config["work_dir"]+"/degeprime_output/Selected_"+config["Output_suffix_name"], p=config["primer_selection_parameters"]
  threads: 1
  message: "Selecting the best primer options after degeprime"
  shell: "python src/Primer_screen.py -d {params.dir} -o {output} {params.p}"

rule primer_screen_fasta:
  input: expand("potential_primers/primers_{suffix}", suffix=config["Output_suffix_name"])
  output: expand("potential_primers/best_primers_pre_selection_{suffix}.fasta", suffix=config["Output_suffix_name"])
  threads: 1
  message: "Printing out the best primer options after degeprime in fasta format"
  shell: "bash src/primer_sequences.sh {input} {output}"

rule mfe_hairpin:
  input: expand("potential_primers/best_primers_pre_selection_{suffix}.fasta", suffix=config["Output_suffix_name"])
  output: expand("potential_primers/hairpin_best_primers_{suffix}", suffix=config["Output_suffix_name"])
  threads: config["threads"]
  message: "Checking hairpin"
  shell: "src/MFEprimer3/mfeprimer hairpin -i {input} -o {output} -c {threads}"

rule hairpin_screen:
  input: i=expand("potential_primers/hairpin_best_primers_{suffix}", suffix=config["Output_suffix_name"]),
         f=expand("potential_primers/best_primers_pre_selection_{suffix}.fasta", suffix=config["Output_suffix_name"])
  output: expand("potential_primers/primers_after_hairpin_screen_{suffix}.fasta", suffix=config["Output_suffix_name"])
  message: "Filtering primers after hairpin calculation"
  threads: 1
  shell: "python src/parse_mfeprimer_hairpin.py -i {input.i} -f {input.f} -o {output}"

rule mfe_dimer:
  input: expand("potential_primers/primers_after_hairpin_screen_{suffix}.fasta", suffix=config["Output_suffix_name"])
  output: expand("potential_primers/dimer_best_primers_{suffix}", suffix=config["Output_suffix_name"])
  message: "Dimer check on pre_selected primers"
  threads: config["threads"]
  shell: "src/MFEprimer3/mfeprimer dimer -i {input} -o {output} -c {threads}"

rule dimer_screen:
  input: i=expand("potential_primers/dimer_best_primers_{suffix}", suffix=config["Output_suffix_name"]),
         f=expand("potential_primers/primers_after_hairpin_screen_{suffix}.fasta", suffix=config["Output_suffix_name"])
  output: t=expand("Results/List_{suffix}.txt", suffix=config["Output_suffix_name"]),
          ta=expand("Results/primers_{suffix}_table", suffix=config["Output_suffix_name"])
  params: min=config["amplicon_size"].split(",")[0], max=config["amplicon_size"].split(",")[1], gene=config["NCBI_gene_name"]
  message: "Printing out the selected primers"
  shell: "python src/parse_mfeprimer_dimer.py -i {input.i} -f {input.f} -o {output.t} -t {output.ta} -m {params.min} -M {params.max} -s {params.gene}"
