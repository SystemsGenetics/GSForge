import json
import inspect

__all__ = [
    "get_by_json_attr",
    "filter_by_json_attr",
    "params_to_json",
    "kwargs_overlap",
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
    """Build a .json string of the string parameters from a given instance of a
    `param.Parameterized` subclass.

    :param parameterized_instance:
    :param skip:
    :return:
    """
    values = {key: value for key, value in parameterized_instance.get_param_values()
              if isinstance(value, str) and value not in skip}
    return json.dumps(values)


def kwargs_overlap(paramaterized, func):
    """Gets the intersection between the parameters of a paramaterized model and
     the keyword arguments of a function, and returns the current intersection
     and their values as a dictionary."""
    key_set = set(inspect.signature(func).parameters.keys()
                  ).intersection(set(paramaterized.param.objects().keys()))
    return {key: getattr(paramaterized, key) for key in key_set}
