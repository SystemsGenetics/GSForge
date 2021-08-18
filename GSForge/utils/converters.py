"""
Converter functions are those that transform data for IO to / from  other programs.
"""

import copy

import pandas as pd

from ..models._GeneSetCollection import GeneSetCollection


# def dict_to_funce_df(coll_dict):
#     modules = [pd.DataFrame({'genes': values, 'module': name})
#                for name, values in coll_dict.items()]
#     return pd.concat(modules, ignore_index=True)
#
#
#
# def parse_specification(specification: dict = None):
#     pass


def gene_set_collection_to_support_dataframe(collection: GeneSetCollection, specification: dict = None) -> pd.DataFrame:
    """

    Parameters
    ----------
    collection : GeneSetCollection
        A ``GeneSetCollection`` object to be converted to a ``pandas.Dataframe``.

    specification : dict
        TODO: Declare required format of specification dictionary.

    Returns
    -------
    pd.DataFrame
        A dataframe with columns ``['Gene', 'Collection']``.

    """
    if specification is None:
        # Prepare the default set combinations.
        pass

    standard_geneset_spec = specification.get('geneset_keys')

    standard_include, standard_exclude = None, None
    if standard_geneset_spec:
        standard_include = standard_geneset_spec.get('include')
        standard_exclude = standard_geneset_spec.get('exclude')

    # Build standard key-based modules.
    # module_dict = copy.deepcopy(gs_coll.as_dict(keys=standard_include, exclude=standard_exclude))
    module_dict = collection.as_dict(keys=standard_include, exclude=standard_exclude)
    standard_keys = tuple(module_dict.keys())

    # Build standard union, intersection and difference modules, and add them
    # to module_dict, which will store all the modules we create in this function.
    module_dict['standard_intersection'] = collection.intersection(keys=standard_include, exclude=standard_exclude)
    module_dict['standard_union'] = collection.union(keys=standard_include, exclude=standard_exclude)

    for primary_key in standard_keys:
        module_dict[f'std_diff__{primary_key}'] = collection.difference(primary_key=primary_key)

    # Build custom intersection modules.
    custom_intersection_specs = specification.get('intersections')
    if custom_intersection_specs:
        for intersection_sepc in custom_intersection_specs:
            module_dict[f"cust_intersect__{intersection_sepc['name']}"] = collection.intersection(
                keys=intersection_sepc.get('include'), exclude=intersection_sepc.get('exclude'))

    # Build union modules.
    custom_union_specs = specification.get('unions')
    if custom_union_specs:
        for union_spec in custom_union_specs:
            module_dict[f"cust_union__{union_spec['name']}"] = collection.union(
                keys=union_spec.get('include'), exclude=union_spec.get('exclude'))

    # Build joint differences.
    joint_differences = specification.get('joint_differences')
    if joint_differences:
        for spec in joint_differences:
            module_dict[f"joint_diff__{spec['name']}"] = collection.joint_difference(
                primary_keys=spec.get('primary_keys'), other_keys=spec.get('other_keys'))

    # Combine all of these into one .csv file for FUNCE.
    module_df = dict_to_funce_df(module_dict)

    return module_df
