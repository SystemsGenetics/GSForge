from sklearn.datasets import make_multilabel_classification
import pandas as pd

from ..models import AnnotatedGEM


# TODO: Add count distribution control.
def generated_random_gem(kind="default", **kwargs):
    options = {
        "default": dict(),
    }

    data, labels = make_multilabel_classification(**options[kind])
    return AnnotatedGEM.from_pandas(pd.DataFrame(data), pd.DataFrame(labels),
                                    name="Generated GEM")
