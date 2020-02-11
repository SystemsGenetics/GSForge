import param

from ._Interface import Interface


class OperationInterface(Interface, param.ParameterizedFunction):
    """
    Abstract class for a GEMOperation.

    Every GEMOperation undergoes some argument parsing, then calls self.process(),
    which must be implemented by implemented classes.
    """

    def process(self):
        """Abstract process."""
        raise NotImplementedError

    def __call__(self, *args, **params):
        super().__init__(*args, **params)
        return self.process()
