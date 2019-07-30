# InteractionScores

Dear reader,

I am glad you're here.

Here's a simple procedure that we developped to predict pairwise cell-to-cell interactions from population level gene expression data.

The method considers each cell-cell interaction as a quantitative variable dependent on the relative expression values of the receptor and the ligand. The method works readily and robustly across different data sets and cell types because the expression values are encoded by their rank. The rank-based comparison generates a global and balanced picture of cell-to-cell interactomes without overemphasizing genes with extremely high expression values

For each pairs of cells, and each interaction in our Receptor Ligand DataBase (or any other RL database of your choice), the script ranks the receptor and ligand  and then normalizes the sum of these two ranks so as to get a score ranging from 0 to 1. 

In a context where several pairs of interacting partners are present, the sets of interactions scores can be compared between pairs of cells  with standard tools of statistical learning (e.g., PCA, clustering) to extract cell-type-specific features of the cell-to-cell interactome 

We propose to use z-scores to identify interactions which are specific to certain pairs of interactors. 

We first used this method to predict potential interactions between hematopoietic progenitors and the bone marrow niche in mice, the paper is now in press in the journal Blood.


Please do not hesitate to contact me if you have any questions regarding the method.

Cheers,


Adrien 
