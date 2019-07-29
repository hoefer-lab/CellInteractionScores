# InteractionScores

Dear reader,

I am glad you're here.

Here's a simple procedure that I developped to infer pairwise interactions between cell populations. it is designed to work with bulk gene expression data of cell populations. The whole idea is that the relative expression of a gene is more informative as to the importance of this gene than the actual expression value.

For each pairs of cells, the script ranks the receptors and ligands and then sums the ranks of receptors and ligands and normalizes this score to a 0:1 range so that the interactions are scored from 0 to 1.

I first used this method to predict potential interactions between hematopoietic progenitors and the bone marrow niche in mice, the paper is now in press in the journal Blood.

Cheers,

Adrien
