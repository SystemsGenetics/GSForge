import json

__all__ = [
    "get_by_json_attr",
    "filter_by_json_attr",
    "params_to_json",
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
