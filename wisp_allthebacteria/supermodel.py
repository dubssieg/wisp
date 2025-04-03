import logging
from pathlib import Path
from dataset import Dataset

from model import XGBoostModel
import copy
from utils import cpu_count, get_current_datetime_string, merge_dicts, SystemStatsLogger
from taxdb import TaxDB
from database import Database

LOG = logging.getLogger(__name__)


class SuperModel:
    def __init__(
        self, model_conf: dict, supermodel_conf: dict, database: Database, taxdb: TaxDB
    ):
        self._model_conf = model_conf
        self._supermodel_conf = supermodel_conf
        self._database = database
        self._taxdb = taxdb
        self._dataset = Dataset(database=self._database, taxdb=self._taxdb)
        self._dt = get_current_datetime_string()
        self._routing = self._get_models_routing()

    def train(self):
        for model_tid in self._routing.keys():
            self._get_trained_model(model_tid)

    def _train_model(self, model_tid: int | None):
        with SystemStatsLogger(interval=30):  # FIXME: main.py
            model = self._get_model(model_tid)
            conf = self._routing[model_tid]["conf"]

            train_batch_count = conf["train_batch_count"]
            if train_batch_count == "max":
                train_batch_count = None
            model_path = self._get_model_paths(model_tid)["model"]

            model.train(
                save_path=model_path,
                train_batch_count=train_batch_count,
                eval_patience=conf["eval_patience"],
                eval_batch_count=conf["eval_batch_count"],
                test_batch_count=conf["test_batch_count"],
                num_boost_round=conf["num_boost_round"],
            )

            model.stop()

    def _get_model(self, model_tid: int | None) -> XGBoostModel:
        workspace_path = self._get_model_paths(model_tid)["workspace"]
        conf = self._routing[model_tid]["conf"]
        rank = self._routing[model_tid]["rank"]

        normalize = conf["normalize"]
        if not normalize:
            normalize = None

        return XGBoostModel(
            rank=rank,
            database=self._database,
            normalize=normalize,
            batch_size=conf["batch_size"],
            taxdb=self._taxdb,
            generator_threads=cpu_count(conf["generator_threads"]),
            sample_balance_factor=conf["sample_balance_factor"],
            batch_balance_factor=conf["batch_balance_factor"],
            min_samples_per_class=conf["min_samples_per_class"],
            max_buffer_total_size=conf["max_buffer_total_size"],
            workspace_path=workspace_path,
            start_generator=False,
        )

    def _get_trained_model(self, model_tid: int | None) -> XGBoostModel:
        model = self._load_model(model_tid)
        if model is None:
            self._train_model(model_tid)
            model = self._load_model(model_tid)
        return model

    def _load_model(self, model_tid) -> XGBoostModel | None:
        model = self._get_model(model_tid)
        model_path = self._get_model_paths(model_tid)["model"]
        try:
            model.load(model_path)
            return model
        except FileNotFoundError:
            return None

    def _get_model_paths(self, model_tid: int | None) -> dict:
        conf = self._routing[model_tid]["conf"]
        rank = self._routing[model_tid]["rank"]
        if conf["force_models_dir"]:
            model_path = Path(conf["force_models_dir"]).resolve()
        else:
            model_path = Path(conf["force_models_dir"]).resolve() / self._dt
        if conf["force_workspaces_dir"]:
            workspace_path = Path(conf["force_workspaces_dir"]).resolve()
        else:
            workspace_path = Path(conf["default_workspaces_dir"]).resolve() / self._dt
        workspace_path = workspace_path / rank / str(model_tid)
        model_path = model_path / rank / str(model_tid)
        return {"model": model_path, "workspace": workspace_path}

    def _get_models_routing(self) -> dict:
        LOG.info("Building model routing...")

        def _models_routing(model_tid: int, step_id: int, parent_filter: dict) -> dict:

            step = self._supermodel_conf[step_id]
            conf = self._step_conf(step["config"])
            rank = step["rank"]

            tax_ids = self._dataset._get_tax_ids_by_rank(
                rank=rank,
                min_samples=conf["min_samples_per_class"],
                parent_filter=parent_filter,
            )

            routing[model_tid] = {
                "rank": rank,
                "parent_filter": parent_filter,
                "conf": conf,
                "classes": list(tax_ids.keys()),
            }

            if step_id + 1 < len(self._supermodel_conf):
                for rank_tid in tax_ids.keys():
                    _models_routing(
                        model_tid=rank_tid,
                        step_id=step_id + 1,
                        parent_filter={rank: rank_tid},
                    )

        routing = {}
        _models_routing(model_tid=None, step_id=0, parent_filter=None)
        LOG.debug("Model routing - DONE")
        return routing

    def _step_conf(self, step_conf: dict) -> dict:
        conf = copy.deepcopy(self._model_conf)
        return merge_dicts(conf, step_conf)

    # DEPRECATED FOR NOW
    def _build_tree(self) -> dict:
        def recursive_build(rank_index: int, parent_filter: dict) -> dict:
            if rank_index >= len(self._supermodel_conf):
                return {}

            step = self._supermodel_conf[rank_index]
            conf = self._step_conf(step["config"])
            rank = step["rank"]

            tax_ids = self._dataset._get_tax_ids_by_rank(
                rank=rank,
                min_samples=conf["min_samples_per_class"],
                parent_filter=parent_filter,
            ).keys()

            return {
                tid: recursive_build(
                    rank_index=rank_index + 1, parent_filter={rank: tid}
                )
                for tid in tax_ids
            }

        return recursive_build(0, {})

    def _flatten_tree(
        self, tree: dict, skip_empty: bool = False, skip_single: bool = False
    ) -> dict:
        flat_dict = {}

        def recursive_flatten(node: dict):
            for parent_id, children in node.items():
                flat_dict[parent_id] = list(children.keys())
                recursive_flatten(children)

        recursive_flatten(tree)

        if skip_empty:
            flat_dict = {k: v for k, v in flat_dict.items() if len(v) > 0}

        if skip_empty:
            flat_dict = {k: v for k, v in flat_dict.items() if len(v) > 1}
        return flat_dict

    def _get_model_conf(self, tax_id: int) -> dict:
        rank = self._taxdb[tax_id]["Rank"]
        conf = self._model_conf
        for step in self._supermodel_conf:
            if step["rank"] == rank:
                conf = self._step_conf(step["config"])
                break
        conf["force_models_dir"]
        return {"rank": rank, "conf": conf}
