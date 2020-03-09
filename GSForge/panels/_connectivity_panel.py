# import param
# import panel as pn
# import holoviews as hv
# from methodtools import lru_cache
#
# import networkx as nx
# import itertools
#
# from ..models import Interface
# from GSForge.utils._panel_utils import generate_help_pane
#
#
# class Connectivity_Panel(Interface):
#     """
#     Display lineament collection membership via a network graph.
#
#     This panel requires a GeneSetCollection as input.
#     """
#     mapping_selector = param.ListSelector(default=None)
#     mapping_index_name = param.String(default="Gene", precedence=-1.0)
#     edge_weight_label = param.String(default=None, precedence=-1.0)
#     update_network = param.Action(lambda self: self.param.trigger('update_network'))
#
#     def __init__(self, *args, **params):
#         super().__init__(*args, **params)
#         if self.lineament_collection is None:
#             raise Warning("Requires a `GeneSetCollection` input to function.")
#         if self.mapping_selector is None:
#             avail_mappings = list(self.lineament_collection.lineaments.keys())
#             self.param["mapping_selector"].objects = avail_mappings
#             self.set_param(mapping_selector=avail_mappings)
#
#     def build_nx_graph(self, selected_mappings=None, weight=None):
#         """Construct the networkx graph object form the selected mappings."""
#         lcoll = self.lineament_collection
#         graph = nx.Graph()
#         gene_mappings = lcoll.as_dict(selected_mappings)
#
#         for key, genes in gene_mappings.items():
#             graph.add_node(key, node_type="Feature")
#             graph.add_nodes_from(genes, node_type="Gene")
#
#             if weight is not None:
#                 subset = self.lineament_collection[key].data.sel({self.mapping_index_name: genes})
#                 weights = subset[weight].values
#                 new_edges = [(key, gene, {"weight": weight})
#                              for gene, weight in zip(genes, weights)]
#             else:
#                 new_edges = itertools.product([key, ], genes)
#
#             graph.add_edges_from(new_edges)
#
#         return graph
#
#         # # Get the union of pairwise-intersections of genes.
#         # # This will effectively be all genes of connectivity greater than one.
#         # pw_gene_union = set.union(*[set(x) for x in lcoll.pairwise_intersection(selected_mappings).values()])
#         #
#         # for coll_name, coll_genes in lcoll.as_dict(selected_mappings).items():
#         #     # Add the lineament name as a node.
#         #     graph.add_node(coll_name, node_type="Feature")
#         #     # Filter the genes, keep only those that should have more than one edge.
#         #     kept_genes = set.intersection(set(coll_genes), pw_gene_union)
#         #     graph.add_nodes_from(kept_genes, node_type="Gene")
#         #
#         #     # Now get the genes to be combined into a mega-node. These are genes that
#         #     # should just have an edge to this collection spoke node.
#         #     bundled_genes = list(set.difference(set(coll_genes), kept_genes))
#         #     bundle_node = f"{coll_name} bundle of {len(bundled_genes)} genes"
#         #     graph.add_node(bundle_node, node_type="Bundled Genes", size=len(bundled_genes))
#         #
#         #     if weight is not None:
#         #         subset = self.lineament_collection[coll_name].data.sel({self.mapping_index_name: list(kept_genes)})
#         #         weights = subset[weight].values
#         #         new_edges = [(coll_name, gene, {"weight": weight}) for gene, weight in zip(kept_genes, weights)]
#         #     else:
#         #         new_edges = list(itertools.product([coll_name, ], kept_genes))
#         #
#         #     new_edges += [(coll_name, bundle_node)]
#         #     graph.add_edges_from(new_edges)
#
#         # return graph
#
#     @lru_cache()
#     def layout_graph(self, selected_mappings, weight=None):
#         graph = self.build_nx_graph(selected_mappings, weight)
#         # Construct a circular layout for the collection source nodes.
#         # mapping_layout = nx.circular_layout(selected_mappings)
#         layout = nx.layout.kamada_kawai_layout(graph, weight=None)
#         return graph, layout
#
#     @param.depends("update_network")
#     def view(self):
#         # Create an immutable tuple of mapping keys so that they can be hashed.
#         hashable_keys = tuple(self.mapping_selector)
#         graph, layout = self.layout_graph(hashable_keys, self.edge_weight_label)
#         return hv.Graph.from_networkx(graph, positions=layout).opts(
#             tools=['hover'], padding=0.2, height=400, width=500,
#             color_index="node_type", node_size=5, edge_line_width=.5,
#             xaxis=None, yaxis=None)
#
#     def panel(self):
#         controls = pn.Param(self.param)
#         main_layout = pn.Column("### GeneSet Connectivity Graph",
#                                 pn.Row(self.view, controls))
#         help_layout = generate_help_pane(self)
#         return pn.Tabs(("UMAP", main_layout), ("Documentation", help_layout))
