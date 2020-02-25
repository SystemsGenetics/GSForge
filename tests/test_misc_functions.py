
"""
Test Features Methods GSForge.
"""

import GSForge as gsf


def test_get_gene_index(random_annotated_gem):

	# Create the control interface and the test interface
	inter = gsf.Interface(gem=random_annotated_gem)
	test_inter = gsf.TestInterface(gem=random_annotated_gem)

	# Get both the index values from control and test interface
	index = inter.get_gene_index()
	test_index = test_inter.test_get_gene_index()

	# Compare index values
	for x in range(len(index)):
		assert index[x] == test_index[x]


def test_get_sample_index(random_annotated_gem):

	# Create the control interface and the test interface
	inter = gsf.Interface(gem=random_annotated_gem)
	test_inter = gsf.TestInterface(gem=random_annotated_gem)

	# Get both the index values from control and test interface
	index = inter.get_sample_index()
	test_index = test_inter.test_get_sample_index()

	# Compare index values
	for x in range(len(index)):
		assert index[x] == test_index[x]


def test_gene_support(random_gene_set):
	geneset = gsf.GeneSet(random_gene_set)
	test_geneset = gsf.TestGeneSet(gene_sets=random_gene_set)

	support = geneset.gene_support()
	test_support = test_geneset.test_gene_support()

	# Compare index values
	assert support == test_support
