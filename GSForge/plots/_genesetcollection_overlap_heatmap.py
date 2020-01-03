import itertools
import holoviews as hv


def genesetcollection_overlap_heatmap(alpha: dict, beta: dict):
    """
    View the overlap (as a heatmap) of two GeneSetCollection dictionaries,
    as provided by gsc.as_dict().

    :param alpha:
    :param beta:
    :return:
    """
    overlap_dim = hv.Dimension("Overlap Count", value_format=lambda x: f"{x}")

    data = [(f"{ak}: {len(av)}", f"{bk}: {len(bv)}", len(set.intersection(set(av), set(bv))))
            for (ak, av), (bk, bv) in itertools.product(alpha.items(), beta.items())]

    heatmap = hv.HeatMap(data, vdims=overlap_dim)
    return heatmap * hv.Labels(heatmap)
