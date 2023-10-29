cd results/3.ASV_tree/
conda activate qiime2-2022.8
i=all_f
# all ASVs
echo "Importing ${i}"
qiime tools import --input-path ${i}.fasta --output-path ${i}.qza --type 'FeatureData[Sequence]'
echo "Aligning ${i}"
qiime alignment mafft --i-sequences ${i}.qza --o-alignment ${i}_aln.qza --p-n-threads 6
echo "Masking ${i}"
qiime alignment mask --i-alignment ${i}_aln.qza --o-masked-alignment ${i}_aln_msk.qza
echo "Building phylo for ${i}"
qiime phylogeny fasttree --i-alignment ${i}_aln_msk.qza --o-tree  ${i}_utree.qza --p-n-threads 6
echo "Rooting ${i}"
qiime phylogeny midpoint-root --i-tree ${i}_utree.qza --o-rooted-tree ${i}_rtree.qza
echo "Exporting ${i}"
qiime tools export --input-path ${i}_rtree.qza --output-path ${i}_tree
mv ${i}_tree/tree.nwk ${i}.nwk
rm -r ${i}_tree
conda deactivate
