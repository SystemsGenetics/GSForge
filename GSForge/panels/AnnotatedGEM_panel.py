import panel as pn
import param
import holoviews as hv


class AnnotatedGEM_Panel():

    annotated_gem = param.Parameter()

    def view(self):
        return self.annotated_gem.__repr__()
