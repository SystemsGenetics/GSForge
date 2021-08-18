import holoviews as hv


# TODO: Re-implement below.
# def plot_label_bars(label_df, max_x=10):
#     bar_list = list()
#
#     for col in label_df.columns:
#
#         if label_df[col].nunique() > max_x:
#             continue
#
#         label_names = label_df.groupby([col]).nunique().index.values
#         label_counts = label_df.groupby([col]).count().max(axis=1).values
#
#         new_bar = hv.Bars(zip(label_names, label_counts),
#                           kdims=[col], vdims=["count"],
#                           group='Label Counts', label=col)
#
#         bar_list.append(new_bar)
#
#     return hv.Layout(bar_list)
