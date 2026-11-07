"""
splitgap
========
A diagnostic library for auditing and preventing Sequence-Neighbour Gap
spatial leakage in biological machine learning cross-validation.
"""

__version__ = "1.0.0"

from .core import (
    compute_nearest_train_distances,
    audit_cv_splits,
    partition_modulo,
    partition_contiguous_blocks,
)
from .recommend import (
    recommend_block_partition,
)

__all__ = [
    "compute_nearest_train_distances",
    "audit_cv_splits",
    "partition_modulo",
    "partition_contiguous_blocks",
    "recommend_block_partition",
]
