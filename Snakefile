import os

# Define main directories
data_dir = "/home/utilisateur/Bureau/single_cell_DEG/data"
output_dir = "/home/utilisateur/Bureau/single_cell_DEG/results"

# Dynamically list all directories in data_dir to create sample names
sample_names = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d))]

# Rule all: Specifies the final outputs to be generated
rule all:
    input:
        # Final outputs from the differential_expression rule
        os.path.join(output_dir, "DE_results.csv"),
        os.path.join(output_dir, "volcano_plot.png"),
        os.path.join(output_dir, "PCA_plot.png")

# Rule to load data, perform quality control, and save the Seurat object
# it return name list of seurat objects as such : 
# seurat_list[[sample_name]] <- seurat_obj
# it takes as input :
#name list of sample directory and sample name as well as output dir to save seurat objects
# call : seurat_objects <- load_and_process_data(data_dirs, output_dir)


rule quality_control:
    input:
        input_dir = data_dir  # Input directory (fixed path)
    output:
        "{output_dir}/{sample}_seurat_object.rds"  # Output file for each sample
    params:
        output_dir = output_dir  # Output directory (fixed path)
    shell:
        """
        Rscript scripts/quality_control_and_preprocessing.R {input.input_dir} {params.output_dir}
        """


# Rule for differential expression analysis
# input : the list of seurat object created in previous step and the output dir for results plots
# call : perform_differential_expression(seurat_list, output_dir)

# Rule for differential expression analysis
rule differential_expression:
    input:
        seurat_objects=expand("{output_dir}/{sample}_seurat_object.rds", output_dir=output_dir, sample=sample_names)
    output:
        de_results=os.path.join(output_dir, "DE_results.csv"),
        volcano_plot=os.path.join(output_dir, "volcano_plot.png"),
        pca_plot=os.path.join(output_dir, "PCA_plot.png")
    params:
        output_dir=output_dir
    shell:
        """
        Rscript --vanilla scripts/differential_expression_analysis.R {input.seurat_objects} {params.output_dir} 
        touch {output.de_results} {output.volcano_plot} {output.pca_plot}
        """

