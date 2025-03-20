import sys
sys.path.append('../..')
from wisp.wisp_light.training.metrics import ConfusionMatrixTracker
from wisp.wisp_light.dataset.refSeqDataset import TAXO_LEVELS

dataset = [
    (
        {'phylum': 'A', 'class': 'C', 'order': 'G', 'family': 'K'},  # Ground truth (gt_taxons)
        {'phylum': 'A', 'class': 'C', 'order': 'G', 'family': 'K'}   # Prediction (correct)
    ),
    (
        {'phylum': 'A', 'class': 'C', 'order': 'H', 'family': 'L'},
        {'phylum': 'B', 'class': 'E', 'order': 'J', 'family': 'P'}  # error phylum
    ),
    (
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'M'},
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'M'}  # correct
    ),
    (
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'N'},
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'N'}  # correct
    ),
    (
        {'phylum': 'B', 'class': 'E', 'order': 'J', 'family': 'O'},
        {'phylum': 'B', 'class': 'E', 'order': 'J', 'family': 'O'}  # correct
    ),
    (
        {'phylum': 'B', 'class': 'E', 'order': 'J', 'family': 'P'},
        {'phylum': 'B', 'class': 'E', 'order': 'J', 'family': 'P'}  # correct
    ),
    (
        {'phylum': 'A', 'class': 'C', 'order': 'H', 'family': 'L'},
        {'phylum': 'A', 'class': 'C', 'order': 'G', 'family': 'K'}  # error order
    ),
    (
        {'phylum': 'A', 'class': 'C', 'order': 'H', 'family': 'L'},
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'M'}  # error class
    ),
    (
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'M'},
        {'phylum': 'A', 'class': 'D', 'order': 'I', 'family': 'N'}  # error family
    )

]
metrics = ConfusionMatrixTracker()

for sample in dataset:
    gt_taxons, pred_taxons = sample
    metrics.update(true_labels=gt_taxons, pred_labels=pred_taxons)

metrics.build_taxonomy_df()

for level in TAXO_LEVELS:
    conf_mat_level = metrics.get_confusion_matrix(level)
    print(f"level {level} :\n", conf_mat_level)

print(metrics.taxonomy_df)
level_base ='order'
level_sep = 'class'
separator_indices = metrics.calculate_separator_indices(level_1=level_sep, level_2=level_base)
print("Indices de séparation :", separator_indices)

conf_mat = metrics.get_confusion_matrix(level_base)
from wisp.wisp_light.visu.plots_tools import plot_conf_mat
plot_conf_mat(conf_mat,  level_base, separator_indices, filename=None)

