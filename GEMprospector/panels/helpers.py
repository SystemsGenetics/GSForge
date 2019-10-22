import panel as pn


def generate_help_pane(self, cutoff=-1.0):
    """
    Generate a panel of the docstrings of the model provided.
    """
    docs = []
    for val in self.param.objects():
        if self.param[val].precedence and self.param[val].precedence >= cutoff:
            docs.append(f"### {self.param[val].label}\n\n" + (self.param[val].doc or ""))
    return pn.pane.Markdown("\n\n".join(docs), width=800)
