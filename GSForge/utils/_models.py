import json
import inspect

# noinspection PyProtectedMember
from ..models._GeneSet import GeneSet
# noinspection PyProtectedMember
from ..models._GeneSetCollection import GeneSetCollection

__all__ = [
    "get_by_json_attr",
    "filter_by_json_attr",
    "params_to_json",
    "kwargs_overlap",
    "merge_collections",
    "merge_and_combine_collections",
]


def get_by_json_attr(xr_dataset, key_pair):
    main_key, json_key = key_pair
    json_str = xr_dataset.attrs[main_key]
    json_attrs = json.loads(json_str)
    return json_attrs.get(json_key)


def filter_by_json_attr(xr_dataset, key_pair, req_value):
    if get_by_json_attr(xr_dataset, key_pair) == req_value:
        return True
    else:
        return False


def params_to_json(parameterized_instance, skip: list = None):
    """
    Builds a .json string of the string parameters from a given instance of a ``param.Parameterized`` subclass.

    :param parameterized_instance:
      An instance of a ``param.Parameterized`` class.

    :param skip:
      List of keys to be skipped.

    :return:
      A json encoded string of the parameters.
    """
    values = {key: value for key, value in parameterized_instance.get_param_values()
              if isinstance(value, str) and value not in skip}
    return json.dumps(values)


def kwargs_overlap(paramaterized, func):
    """
    Gets the intersection between the parameters of a ``param.Parameterized`` class and
    the keyword arguments of a function, and returns the current intersection
    and their values as a dictionary.
    """
    key_set = set(inspect.signature(func).parameters.keys()
                  ).intersection(set(paramaterized.param.objects().keys()))
    return {key: getattr(paramaterized, key) for key in key_set}


def merge_collections(*collections, **params):
    """
    Merge the GeneSetCollection objects provided into a single GeneSetCollection.

    This merges all of the collection gene_set dictionaries, and returns them
    within a single GeneSetCollection object.

    :param collections:
        GeneSetCollection objects, these should be named (via the name parameter)
        for the output dictionary to be interpretable.

    :param params:
        Keyword parameters to be passed to the new GeneSetCollection.

    :returns:
        A new GeneSetCollection.
    """
    new_gene_set_dict = dict()
    for coll in collections:
        new_gene_set_dict = {**coll.gene_sets, **new_gene_set_dict}
    return GeneSetCollection(gene_sets=new_gene_set_dict, **params)


def merge_and_combine_collections(*collections, **params):
    """
    Merge the GeneSets within the given GeneSetCollection objects provided and then
    combine them into a a single GeneSetCollection.

    :param collections:
        GeneSetCollection objects, these should be named (via the name parameter)
        for the output dictionary to be interpretable.

    :param params:
        Keyword parameters to be passed to the new GeneSetCollection.

    :return:
        A new GeneSetCollection.
    """
    new_gene_set_dict = dict()
    for coll in collections:
        gene_sets = list(coll.gene_sets.values())
        new_gene_set_dict[coll.name] = GeneSet.from_GeneSets(*gene_sets)
    return GeneSetCollection(gene_sets=new_gene_set_dict, **params)
