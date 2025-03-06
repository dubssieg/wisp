from collections import defaultdict
from pathlib import Path
import networkx as nx
from tqdm.auto import tqdm


from api import API

RANKS = {
    "biotype",
    "clade",
    "class",
    "family",
    "genus",
    "kingdom",
    "no rank",
    "order",
    "pathogroup",
    "phylum",
    "serogroup",
    "serotype",
    "species",
    "species group",
    "species subgroup",
    "strain",
    "subclass",
    "subfamily",
    "subgenus",
    "suborder",
    "subspecies",
    "superkingdom",
    "tribe",
}


class TaxoGraph:
    def __init__(self, api: API):
        self._api = api
        self._graph = nx.DiGraph()
        self._db_tax_ids = list()

    def build_graph_from_db(self, db_path: str | Path):
        """Build the taxonomy graph from the cache data."""
        db_path = Path(db_path)
        # tax_id list from db_path
        self._db_tax_ids = [
            int(dir.name)
            for dir in db_path.iterdir()
            if dir.is_dir() and dir.name.isdigit()
        ]

        # build graph
        for tax_id in tqdm(self._db_tax_ids):
            tdata = self._api[tax_id]
            self._graph.add_node(
                tax_id, name=tdata["ScientificName"], rank=tdata["Rank"]
            )

            # get full lineage
            lineage = tdata.get("LineageEx", [])
            for i in range(len(lineage) - 1):
                parent = lineage[i]
                child = lineage[i + 1]

                # parent node
                if int(parent["TaxId"]) not in self._graph:
                    parent_data = self._api[int(parent["TaxId"])]
                    self._graph.add_node(
                        int(parent["TaxId"]),
                        name=parent_data["ScientificName"],
                        rank=parent_data["Rank"],
                    )

                # child node
                if int(child["TaxId"]) not in self._graph:
                    child_data = self._api[int(child["TaxId"])]
                    self._graph.add_node(
                        int(child["TaxId"]),
                        name=child_data["ScientificName"],
                        rank=child_data["Rank"],
                    )

                # parent - child node
                self._graph.add_edge(int(parent["TaxId"]), int(child["TaxId"]))

    def get_ranks_by_tax_id(self, rank: str) -> dict:
        rank_mapping = defaultdict()
        for tax_id in self._db_tax_ids:
            lineage_ex = self._api[tax_id].get("LineageEx", [])
            for entry in lineage_ex:
                if entry["Rank"] == rank:
                    rank_mapping[tax_id] = entry["TaxId"]
                    break

        return rank_mapping

    # def get_tax_ids_by_rank(self, rank: str) -> dict:
    #     """Return all tax_ids classified by the given rank that were in db_path."""
    #     # find nodes with the specified rank
    #     rank_tax_ids = [
    #         node
    #         for node, attrs in self._graph.nodes(data=True)
    #         if attrs["rank"] == rank
    #     ]

    #     # Get all descendants of these nodes
    #     tax_ids_by_rank = dict()
    #     for tax_id in rank_tax_ids:
    #         descendants = set(nx.descendants(self._graph, tax_id))
    #         descendants = descendants.intersection(self._db_tax_ids)
    #         tax_ids_by_rank[tax_id] = list(descendants)
    #     return tax_ids_by_rank

    # def get_ancestors(self, tax_id: int) -> list:
    #     """Return all ancestors of a given tax_id."""
    #     return list(nx.ancestors(self._graph, tax_id))
