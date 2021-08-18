import numpy as np
from scipy.stats import fisher_exact


def count_set_terms(query_genes, term_df, min_term_count=2):
    gene_annots = fdf.groupby('Gene')['Term'].unique()
    gene_overlap = np.intersect1d(gene_annots.index, query_genes)
    term_counts = term_df[term_df['Gene'].isin(gene_overlap)]['Term'].value_counts()
    return term_counts[term_counts > min_term_count]


def get_genes_by_term(query_term, term_df):
    return np.unique(np.concatenate(term_df.groupby('Term')['Gene'].unique()[query_term].values))


def enrichment_score(query_genes, term_genes, background_genes):
    """Calculalte an enrichment score of term_genes in query_genes with a given background.
    
        |                | Has term     | Not Has term |
        | ---------------| -------------| ------------ |
        | In module      |   n11        |    n12       |
        | Not in module  |   n21        |    n22       |
    """
    n11 = np.intersect1d(query_genes, term_genes).size
    n12 = query_genes.size - n11
    n21 = term_genes.size - n11
    n22 = background_genes.size - (n11 + n12 + n21)
    contmatrix = np.array([[n11, n12], [n21, n22]])
    odds_ratio, p_values = stats.fisher_exact(contmatrix, alternative='greater')
    return odds_ratio, p_values


def ks_enrichment_score(gene_index, gene_scores, term_genes, exp=1):
    sorted_genes = gene_index[np.argsort(gene_scores)[::-1]]  # Sort high to low.
    hits = np.isin(gene_index, term_genes)
    miss = np.invert(hits)
    hits = hits.astype(int)
    miss = miss.astype(int)
    n_genes = gene_index.shape[0]
    n_terms = term_genes.shape[0]
    hit_sum = np.sum(gene_scores ** exp * hits)
    hits_scores = np.cumsum(gene_scores ** exp * hits / hit_sum)
    miss_scores = np.cumsum(miss / (n_genes - np.sum(hits)))
    scores = hits_scores - miss_scores
    stat, pval = stats.ks_2samp(hits_scores, miss_scores, 'greater')
    return scores, stat, pval


# def ks_score_matrix(measure, background_index, hit_index) -> np.ndarray:
#     """Performs a Kolmogorov-Smirnov type test for enrichment calculation based on a ranking measure.
    
#     Parameters
#     ----------
#     measure : np.ndarray
#         A pairwise distance matrix.

#     background_index : np.ndarray
#         The background index, should work for row and column labels of the measure parameter.
        
#     hit_index : np.ndarray
#         The index values that are to be checked for enrichment.
    
#     Returns
#     -------
#     scores : np.ndarray
#         A matrix where each row represents a index value.
#         The row values are cummulative sum of the sorted measures.
#     """
#     n = background_index.shape[0]
#     nt = hit_index.shape[0]
    
#     sorted_idx = np.argsort(measure, axis=1)[::-1]  # Sort high to low.
#     term_pos = np.isin(background_index, hit_index)
    
#     # Use the sort indexes to transform the term position 
#     # boolean array into the 'hit' matrix.
#     hits = term_pos[sorted_idx]
#     miss = np.invert(hits)
    
#     # Sum the hit scores for each gene.
#     hit_sum = np.sum(measure * hits, axis=1)
    
#     hits_scores = np.cumsum(measure * hits / hit_sum, axis=1)
#     miss_scores = np.cumsum(miss / (n - nt), axis=1)
#     scores = hits_scores - miss_scores
#     return scores

