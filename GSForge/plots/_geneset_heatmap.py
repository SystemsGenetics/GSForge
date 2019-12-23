import holoviews as hv

import itertools


def geneset_overlap_heatmap(lineament_collection, keys=None, mode="overlap"):
    lineament_dict = lineament_collection.as_dict(keys)

    if mode == "overlap":
        overlap_dim = hv.Dimension("Overlap Count", value_format=lambda x: f"{x}")
        data = [(ak, bk, len(set.intersection(set(av), set(bv))))
                for (ak, av), (bk, bv) in itertools.permutations(lineament_dict.items(), 2)
                if ak != bk]

    elif mode == "percent":
        overlap_dim = hv.Dimension("Overlap Percent", value_format=lambda x: f"{x:.0%}")

        zero_filtered_dict = {k: v for k, v in lineament_dict.items()
                              if len(v) > 0}
        data = [(ak, bk, len(set.intersection(set(av), set(bv))) / len(set(av)))
                for (ak, av), (bk, bv) in itertools.permutations(zero_filtered_dict.items(), 2)
                if ak != bk]

    else:
        raise ValueError(f"{mode} is not a valid mode. Select from 'overlap' or 'percent'.)")

    heatmap = hv.HeatMap(data, vdims=overlap_dim)  # .options(xrotation=45, width=450, height=450, labelled=[],
    #         colorbar=True)
    return heatmap * hv.Labels(heatmap)
