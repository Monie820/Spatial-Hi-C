
## Spatial Hi-C preprocessing workflow
We used [Higashi](https://github.com/ma-compbio/Higashi) to impute single spot contacts matrices, and performed single spot A/B compartment calling. In details,
1. Cropped spots within the tissue region [seruat_HiC_crop.R](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/seruat_HiC_crop.R).
2. Fastq are processed, mapped, convert to bedpe file and prepare higashi file using [base_pipeline_HPC_MboI_mouse_1_50.sh](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/base_pipeline_HPC_MboI_mouse_1_50.sh).
3. High quality spots are selected using the notebook [extract_barcode_based_on_knees.R](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/extract_barcode_based_on_knees.R).
4. Contacts information from individual cells are extracted using scripts [stat_cropped_from_bedpe.py](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/stat_cropped_from_bedpe.py).
5. Single spot A/B compartment score is calculated using Higashi by [Higashi_h_500kb.ipynb](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/Higashi_h_500kb.ipynb). This allows us to perform embedding and spot annotation.

## Spatial transcriptome preprocessing workflow
Fastq are processed, mapped, and convert to gene expression matrix file using [stpipeline.sh](https://github.com/Monie520/Spatial-Hi-C/blob/main/1.preprocess/stpipeline.sh).
