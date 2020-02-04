import panel as pn
import param
import holoviews as hv


class InterfacePanel(param.Parameterized):

    interface = param.Parameter()

    def view(self):
        count_table = hv.Table(self.interface.x_count_data)
        return count_table

    def panel(self):
        pass