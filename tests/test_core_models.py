"""
Test the 'core' data models, i.e. those that have ``.data`` parameter
with an ``xarray.Dataset`` object.
"""
import GSForge as gsf


def test_annotated_gem_creation(random_annotated_gem):
    """ Test if an AnnotatedGEM can be constructed. """
    assert random_annotated_gem


def test_gene_set_creation(random_gene_set):
    """Test if a GeneSet can be constructed."""
    assert random_gene_set


def test_gene_set_collection_creation(random_gene_set_collection):
    print(random_gene_set_collection.__repr__())
    as_dict = random_gene_set_collection.as_dict()

    print(".as_dict() output:")
    for key, values in as_dict.items():
        print(key, values.shape)

    set_operations = dict(
        intersection=random_gene_set_collection.intersection(),
        union=random_gene_set_collection.union(),
        difference=random_gene_set_collection.difference())

    print("set operation shapes:")
    for key, values in set_operations.items():
        print(key, values.shape)


def test_interface(random_gene_set_collection):
    interface = gsf.Interface(random_gene_set_collection)
    interface.get_gene_index()


