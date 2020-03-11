import copy
import functools
import json
import os
from pathlib import Path
from textwrap import dedent
from typing import Union, AnyStr, IO

import numpy as np
import pandas as pd
import param
import xarray as xr

from ._GeneSet import GeneSet

class TestGeneSet(GeneSet):

    def test_gene_support(self) -> np.ndarray:
        """
        Returns the list of genes 'supported in this GeneSet.

        The value that this return is (by default) controlled by the self.support_index_name
        parameter.

        Returns
        -------
            A numpy array of the genes 'supported' by this GeneSet.
        """
        if self.support_index_name in self.data:
            supported_genes = self.data[self.gene_index_name].sel(
                {self.gene_index_name: (self.data[self.support_index_name] == True)}).values.copy()
            return self.data[self.gene_index_name].sel({self.gene_index_name: supported_genes}).values.copy()
        else:
            Warning("Warning, no boolean support array found. Returning the complete "
                    "set of genes in this GeneSet.")
            return self.data[self.gene_index_name].values.copy()

    @property
    def test_gene_index(self) -> xr.DataArray:
        """
        Returns the entire gene index of this GeneSet object as an ``xarray.DataArray``.

        The variable or coordinate that this returns is controlled by the `gene_index_name`
        parameter.

        Returns
        -------
        xr.DataArray : A copy of the entire gene index of this GeneSet as an ``xarray.DataArray``.
        """
        return self.data[self.gene_index_name].copy()

    def test_get_n_top_genes(self, score_variable: str = None, mode: str = "absolute_largest",
                        within_support: bool = True, **mode_kwargs) -> np.ndarray:
        """
        Get the top n genes as an array, selected via `mode` by using `score_variable` from the GeneSet.data variables.

        Parameters
        ----------
        score_variable :
            The name of the score variable to use. If this is not given some common score variable names are
            used, and the first found is used.

        mode : str
            The mode by which to calculate how genes are ranked.

            **Available modes**:

            `absolute_largest`
                Returns the n absolute largest values per the `score_variable` variable.
                Useful for selecting the absolute largest log-fold-change, regardless of the sign.

        within_support : bool
            Whether to limit the genes returned to those within the "support" of this GeneSet. Defaults to True.

        mode_kwargs : dict
            Remaining keyword arguments are passed to the 'mode' function.

        Returns
        -------
            np.ndarray : An array of n, top-ranked, genes.
        """

        # TODO: Comment mode functions and their arguments clearly. Include how to write such a function.
        # TODO: Consider moving these functions to utilities?
        def _abs_largest(scores: xr.DataArray, **kwargs):
            nn = kwargs["n"]
            top_idx = np.argsort(np.abs(scores.values))[::-1][:nn]
            return scores.isel({self.gene_index_name: top_idx})[self.gene_index_name].values.copy()

        def _above_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: scores.values > threshold})[self.gene_index_name].values.copy()

        def _below_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: scores.values < threshold})[self.gene_index_name].values.copy()

        def _above_absolute_threshold(scores: xr.DataArray, **kwargs):
            threshold = kwargs["threshold"]
            return scores.isel({self.gene_index_name: np.abs(scores.values) > threshold})[
                self.gene_index_name].values.copy()

        modes = {
            "absolute_largest": _abs_largest,
            "above_threshold": _above_threshold,
            "below_threshold": _below_threshold,
            "above_absolute_threshold": _above_absolute_threshold,
        }

        def _scan_for_score_variable(target=score_variable):
            score_variables = [
                "score",
                "logFC",  # from EdgeR.
                "log2FoldChange",  # from DESeq2.
                "pvalue",  # from DESeq2.
            ]
            if target is not None:  # An explicit variable has been given.
                return target
            else:
                for name in score_variables:
                    for var_name in list(self.data.variables.keys()):
                        if name == var_name:
                            print(f"Automatically found a score variable, {name}. Ensure this the desired variable.")
                            return name

        score_name = _scan_for_score_variable()

        if score_name is None:
            raise ValueError("A score variable could not be automatically identified, please specify one.")

        data = self.data.sel({self.gene_index_name: self.gene_support()})[score_name] \
            if within_support \
            else self.data[score_name]

        # Assign the selected mode function.
        selected_mode = modes.get(mode)
        if selected_mode is None:
            raise ValueError(f"Given mode selection: {mode} is not a valid mode. Choose from:\n{modes.keys()}")

        return selected_mode(data, **mode_kwargs)