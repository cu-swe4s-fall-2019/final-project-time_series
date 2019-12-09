test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

UP=100
DOWN=100
INPUT=test.tsv
BED=gene_location_test.bed

run good_params python promoter_Finder.py \
	--input_genes $INPUT \
	--bedfile $BED \
	--upstream_distance $UP \
	--downstream_distance $DOWN \
	--fasta test.fa \
	--motif_file HOCOMOCOv11_core_MOUSE_mono_meme_format.meme
assert_exit_code 0

# Bad parameters should throw an error
run bad_params python promoter_Finder.py \
	--input_genes 100 \
	--bedfile $BED \
	--upstream_distance $UP \
	--downstream_distance $DOWN \
	--fasta test.fa
assert_exit_code 1
