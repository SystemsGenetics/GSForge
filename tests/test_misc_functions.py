
"""
Test Features Methods GSForge.
"""

import GSForge as gsf
import numpy as np


def test_get_gene_index(random_annotated_gem):

	# Create the control interface and the test interface
	inter = gsf.Interface(gem=random_annotated_gem)
	test_inter = gsf.TestInterface(gem=random_annotated_gem)

	# Get both the index values from control and test interface
	index = inter.get_gene_index()
	test_index = test_inter.test_get_gene_index()

	# Compare index values
	assert np.array_equal(index, test_index)


def test_get_sample_index(random_annotated_gem):

	# Create the control interface and the test interface
	inter = gsf.Interface(gem=random_annotated_gem)
	test_inter = gsf.TestInterface(gem=random_annotated_gem)

	# Get both the index values from control and test interface
	index = inter.get_sample_index()
	test_index = test_inter.test_get_sample_index()

	# Compare index values
	assert np.array_equal(index, test_index)


def test_y_annotation_data(random_annotated_gem):
	# Create the control interface and the test interface
	inter = gsf.Interface(gem=random_annotated_gem)
	test_inter = gsf.TestInterface(gem=random_annotated_gem)

	y_data = inter.y_annotation_data
	test_y_data = test_inter.test_y_annotation_data

	assert np.array_equal(y_data, test_y_data)


def test_gene_support(random_gene_set):
	df = random_gene_set.to_dataframe(only_supported=False)
	test_geneset = gsf.TestGeneSet(df)

	support = random_gene_set.gene_support()
	print(support)
	test_support = test_geneset.test_gene_support()

	# Compare index values
	assert np.array_equal(support, test_support)


def test_gene_index(random_gene_set):
	df = random_gene_set.to_dataframe(only_supported=False)
	test_geneset = gsf.TestGeneSet(df)

	index = random_gene_set.gene_index
	test_index = test_geneset.test_gene_index

	assert np.array_equal(index.values, test_index.values)


def test_get_n_top_genes(random_gene_set):
	df = random_gene_set.to_dataframe(only_supported=False)
	test_geneset = gsf.TestGeneSet(df)

	# "threshold"
	# score variable "scores"

	supports = [True, False]
	modes = ["absolute_largest",
            "above_threshold",
            "below_threshold",
            "above_absolute_threshold"]
	threehold = .5
	score_variable = "scores"

	for mode in modes:
		for support in supports:
			if mode == "absolute_largest":
				for n in range(0, 100, 1):
					genes = random_gene_set.get_n_top_genes(score_variable=score_variable, mode=mode, within_support=support, n=n)
					test_genes = test_geneset.get_n_top_genes(score_variable=score_variable, mode=mode, within_support=support, n=n)
					assert np.array_equal(genes, test_genes)
			else:
				genes = random_gene_set.get_n_top_genes(score_variable=score_variable, mode=mode, within_support=support, threshold=threehold)
				test_genes = test_geneset.get_n_top_genes(score_variable=score_variable, mode=mode, within_support=support, threshold=threehold)
				assert np.array_equal(genes, test_genes)
