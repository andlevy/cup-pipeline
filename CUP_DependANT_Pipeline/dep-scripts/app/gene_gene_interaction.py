import igraph as ig
import pandas as pd

from app.tools.file_paths import StaticPaths
from app.tools.logger import create_logger
from app.tools.exceptions import BackendError

logger = create_logger()


class GGI:

    def __init__(self, scores: pd.DataFrame):
        self.scores: pd.DataFrame = scores
        self.ggi: ig.Graph = GGI._get_ggi(scores)

    @staticmethod
    def __get_edges():
        logger.debug(f"Reading reactome file from {StaticPaths.reactome}")
        try:
            edges = pd.read_csv(StaticPaths.reactome)
            required_edges = ["gene1", "gene2"]
            if not all(param in edges.columns for param in required_edges):
                raise BackendError(f"Reactome file {StaticPaths.reactome} not in correct format - "
                                   f"expected columns [gene1, gene2], received columns {edges.columns}")
            return edges
        except Exception as e:
            logger.exception("Failed to read reactome file")
            raise BackendError(f"Failed to read reactome file {StaticPaths.reactome} using pandas") from e

    @staticmethod
    def _get_ggi(scores: pd.DataFrame) -> ig.Graph:
        logger.debug(f"GGI._get_ggi(scores ~ {scores.shape})")
        edges = GGI.__get_edges()

        logger.info("Reducing reactome edges to genes that are positively expressed and weighting edges.")
        vertices = set([k for k, v in scores.items()
                        if v > 0])
        valid_edges = edges.loc[(edges["gene1"].isin(vertices)) & (edges["gene2"].isin(vertices))].copy()

        logger.info("Creating unique ids for nodes in the network")
        ids = ig.UniqueIdGenerator()
        edge_list = list(zip(valid_edges["gene1"], valid_edges["gene2"]))
        edge_ids = [(ids[x[0]], ids[x[1]]) for x in edge_list]
        weights = valid_edges["gene1"].map(scores).to_list()

        logger.info("Creating a weighted directed network from mutation and expression information")
        ggi = ig.Graph(edge_ids,
                       vertex_attrs=dict(name=ids.values()),
                       edge_attrs=dict(weight=weights),
                       directed=True)
        return ggi

    def get_features(self) -> pd.DataFrame:
        logger.debug("GGI.get_features()")
        weights = [i['weight'] for i in self.ggi.es]
        genes = [i['name'] for i in self.ggi.vs]

        betweenness = self.ggi.betweenness(weights=weights)
        constraint = self.ggi.constraint(weights=weights)
        coreness = self.ggi.coreness(None)
        degree = self.ggi.degree()
        eccentricity = self.ggi.eccentricity()
        eigenvector_centrality = self.ggi.eigenvector_centrality(weights=weights)
        harmonic_centrality = self.ggi.harmonic_centrality()
        hub_score = self.ggi.hub_score(weights=weights)
        neighbourhood_size1 = self.ggi.neighborhood_size(order=1)
        neighbourhood_size2 = self.ggi.neighborhood_size(order=2)
        neighbourhood_size6 = self.ggi.neighborhood_size(order=6)
        pagerank = self.ggi.pagerank(weights=weights)
        transitivity = self.ggi.transitivity_local_undirected(mode="zero")

        index = ['betweenness',
                 'constraint',
                 'closeness',
                 'coreness',
                 "degree",
                 'eccentricity',
                 'eigenvector_centrality',
                 'harmonic_centrality'
                 'hub_score',
                 'neighbourhood_size1',
                 'neighbourhood_size2',
                 'neighbourhood_size6',
                 'pagerank',
                 'transitivity']

        df = pd.DataFrame([betweenness,
                           constraint,
                           coreness,
                           degree,
                           eccentricity,
                           eigenvector_centrality,
                           harmonic_centrality,
                           hub_score,
                           neighbourhood_size1,
                           neighbourhood_size2,
                           neighbourhood_size6,
                           pagerank,
                           transitivity], columns=genes, index=index)
        df = df.apply(lambda row: row / row.sum(), axis=1)
        features = df.T
        return features
