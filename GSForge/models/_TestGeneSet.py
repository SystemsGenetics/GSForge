from __future__ import annotations

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
