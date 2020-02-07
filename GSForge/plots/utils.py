import inspect

default_arguments = dict(
    lfc_var=[
        "logFC",  # from EdgeR.
        "log2FoldChange",  # from DESeq2.
    ],
    mean_var=["baseMean"],
    pval_var=["pvalue"]
)


def infer_default_kwargs(function):
    def new_kwarged_function(data, *args, **kwargs):
        argument_spec = inspect.getfullargspec(function)._asdict()

        data_variables = set(data.variables.keys())
        function_arg_keys = set(inspect.signature(function).parameters.keys())
        arg_key_overlap = set.intersection(function_arg_keys, set(default_arguments.keys()))

        for arg in arg_key_overlap:
            kwargs[arg] = [v for v in set.intersection(data_variables, default_arguments[arg])][0]

        kwargs["data"] = data
        if args and argument_spec.get("args") is not None:
            for arg_name, arg in zip(argument_spec["args"][1:], args[1:]):
                kwargs[arg_name] = arg
        return function(**kwargs)

    return new_kwarged_function
