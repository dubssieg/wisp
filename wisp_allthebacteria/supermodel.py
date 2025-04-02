from dataset import Dataset

# from model import XGBoostModel
import copy
from utils import merge_dicts

"""
kingdom -> phylum -> class -> order -> family

step 1: kingdom -> 4 classes : k1, k2, k3, k4

step 2-k1 : phylum, uniquement les phylum filtrés par k1 -> n classes k1-p1, k1-p2, etc.
step 2-k2 : phylum, uniquement les phylum filtrés par k2 -> n classes k2-p1, k2-p2, etc.
...

step 3 : class, uniquement les classes filtrées par k1-p1 -> n classes k1-p1-c1
...

"""


class SuperModel:
    def __init__(self, dataset: Dataset, conf: dict):
        self._dataset = dataset
        self._conf = conf
        self._steps = []

    def _step_conf(self, step_conf: dict) -> dict:
        conf = copy.deepcopy(self._conf)
        return merge_dicts(conf, step_conf)

    def _build_tree(self, name: str) -> dict:
        def recursive_build(rank_index: int, parent_filter: dict) -> dict:
            if rank_index >= len(self._conf["supermodels"][name]):
                return {}

            step = self._conf["supermodels"][name][rank_index]
            conf = self._step_conf(step["config"])
            rank = step["rank"]

            tax_ids = self._dataset._get_tax_ids_by_rank(
                rank=rank,
                min_samples=conf["model"]["min_samples_per_class"],
                parent_filter=parent_filter,
            ).keys()

            return {
                tid: recursive_build(
                    rank_index=rank_index + 1, parent_filter={rank: tid}
                )
                for tid in tax_ids
            }

        return recursive_build(0, {})


# {
#     "kingdom": {
#         123: {
#             "phylum": {
#                 12301: {
#                     "class": {
#                         1230101: {},
#                         1230102: {},
#                         1230103: {},
#                         1230104: {},
#                         1230105: {},
#                     }
#                 },
#                 12302: {},
#             },
#             12303: {},
#             12304: {},
#         },
#         345: {},
#         789: {},
#     }
# }
