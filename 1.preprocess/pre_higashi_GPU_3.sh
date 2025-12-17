# source /home/server/miniconda3/bin/activate base
. /date/scripts/Pipelines/temp_bashrc

pixel=$3
sample=$2
work_path=$1
script_path=/date/guomaoni/CZP/spatial_hic/scripts_$pixel
bedpe_path=$work_path/5_bedpe
crop_path=$work_path/6_crop
higashi_path=$work_path/7_higashi_input
higashi_out=$work_path/8_higashi_out

if [ -d $work_path ] && [ ! -z $sample ]; then
	# conda activate higashi
	mkdir $higashi_path
	python $script_path/extract_higashi_input.py \
	--input_stat=$crop_path/stat_cropped.csv \
	--all=$bedpe_path/data_all_pixel.txt \
	--out=$higashi_path/data.txt \
	--fig=$higashi_path/heatmap_higashi.png
	
	python $script_path/create_json_and_pickle_GPU.py \
	--projectName=$sample \
	--chromNum=19 \
	--higashiInput=$higashi_path \
	--tempDir=$higashi_out \
	--chromsizeDir=/date/guomaoni/CZP/spatial_hic_rna/HiC/GRCm38.chrom.sizes.txt \
	--cytoDir=/date/guomaoni/CZP/spatial_hic_rna/HiC/GRCm38_cytoband.fromUCSC.txt \
	--pickleDir=$higashi_path/label_info.pickle \
	--jsonDir=$higashi_path/${sample}.json \
	--input_stat=$crop_path/stat_cropped.csv
else
	echo "Please give the mother_path(of the project name) and the project name in order!
	Just like:bash xxx.sh [/media/xxx] [Spatial-E13]
	Or contact with genger@wechat:13958598285"
fi
