from __future__ import annotations

import numpy as np

from ._Interface import Interface


# Defining a test class for the interface class
class TestInterface(Interface):

    def test_get_gene_index(self, count_variable=None) -> np.array:
        """
        Get the currently selected gene index as a numpy array. This function tests if 
        the .copy should be appended to masked counts. Old way is appended.

        :param count_variable:
            Optional, for selecting alternate count arrays. The variable to be retrieved.

        :return:
            A numpy array of the currently selected genes.
        """
        gene_set_combinations = {
            "union": lambda sel_gs: self.gene_set_collection.union(sel_gs),
            "intersection": lambda sel_gs: self.gene_set_collection.intersection(sel_gs),
            "complete": lambda sel_gs: self.gem.gene_index,
        }

        mask_modes = {
            "complete": lambda counts: counts.fillna(0),
            "masked": lambda counts: counts.where(counts > 0.0),
            "dropped": lambda counts: counts.where(counts > 0.0).dropna(dim=self.gem.gene_index_name),
        }

        # Ensure the correct count array is selected.
        count_variable = self.count_variable if count_variable is None else count_variable

        # If an explicit list of genes is provided, use those and move on to masking.
        if self.selected_genes is not None:
            support = self.selected_genes

        # Check if a collection was provided, if not: return the entire gene index for masking.
        elif self.selected_gene_sets == [None] or not self.gene_set_collection:
            support = self.gem.gene_index

        # Otherwise, use  some combination of GeneSet supports should be used.
        else:
            support = gene_set_combinations[self.gene_set_mode](self.selected_gene_sets)

        # Now the counts can be selected based on the gene support.
        counts = self.gem.data.sel({self.gem.gene_index_name: support})[count_variable]

        # Mask the counts.
        masked_counts = mask_modes[self.count_mask](counts)

        # Get the gene index from this mask.
        # A .copy() of the array values should be used to prevent strange index behavior.
        genes = masked_counts[self.gem.gene_index_name].values.copy()
        return genes

    # function that needs testing
    def test_get_sample_index(self) -> np.ndarray:
        """
        Get the currently selected sample index as a numpy array. This function tests if 
        the .copy should be appended to the return value. Old way is appended.

        Returns
        -------
        np.ndarray : An array of the currently selected samples.
        """

        if self.sample_subset is not None:
            subset = self.gem.data.sel({self.gem.sample_index_name: self.sample_subset})
        else:
            subset = self.gem.data

        if self.annotation_mask == "complete":
            pass
        elif self.annotation_mask == "dropped":
            if self.annotation_variables is not None:
                subset = subset[self.annotation_variables].dropna(dim=self.gem.sample_index_name)

            subset = subset.dropna(dim=self.gem.sample_index_name)

        return subset[self.gem.sample_index_name].values.copy()
