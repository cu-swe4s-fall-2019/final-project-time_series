test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

# Clean up files from previous tests if we have them
rm -f C0_promoters.bed
rm -f C1_promoters.bed
rm -f C2_promoters.bed
rm -f C3_promoters.bed
rm -f C4_promoters.bed
rm -f C5_promoters.bed
rm -f C6_promoters.bed

# In order to streamline testing, I'm using some small test
# files in the same format as our real data. This reduces
# testing time while keeping the functional tests robust

run good_params python main.py \
	--input_genes './ref/test.tsv' \
	--bedfile './ref/test_refseq.bed' \
	--genome './ref/test.fa' \
	-m HOCOMOCOv11_core_MOUSE_mono_meme_format.meme
assert_exit_code 0

# we should have generated a directory called clust_out
run file_checker ls
assert_in_stdout clust_out

# and in clust_out we should have a few files
run clust_checker ls clust_out/
assert_in_stdout Clusters_Objects.tsv
assert_in_stdout Clusters_profiles.pdf
assert_in_stdout Summary.tsv

# We also should have generated a couple of promoter files with specific line counts
# We'll check one here to make sure it's the proper size

run file_size_checker wc -l C1_promoters.bed
assert_in_stdout 20

# Within that file, we should have generated a specific FASTA sequence,
# so let's check that output
run fasta_checker cat C1_promoters.bed
assert_in_stdout chr1:1-200
assert_in_stdout gggaagccccc

# Check to make sure that the plots got created as well
run plot_checker ls out/
assert_in_stdout cluster.png
assert_in_stdout gene_heatmap.png
assert_in_stdout gene_trajectory.png

# These plots will have to be manually evaluated by eye to make sure
# they are correct
