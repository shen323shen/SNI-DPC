"""Public estimator interface for SNI-DPC."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._core import mng_dpc


class SNIDPC:
    """Automatic density-peaks clustering with screened natural neighbors.

    SNI-DPC discovers the natural-neighbor order and output cluster count from
    the data. The estimator does not accept the true cluster count or a
    dataset-specific distance cutoff.
    """

    def fit(self, X: ArrayLike) -> "SNIDPC":
        data = _validated_data(X)
        labels, diagnostics = mng_dpc(data)

        self.labels_: NDArray[np.int_] = np.asarray(labels, dtype=int)
        self.diagnostics_: dict[str, Any] = diagnostics
        self.n_clusters_: int = int(diagnostics["num_clusters"])
        self.natural_neighbor_order_: int = int(
            diagnostics["original_info"]["r"]
        )
        self.low_confidence_indices_: NDArray[np.int_] = np.asarray(
            diagnostics["low_indices"],
            dtype=int,
        )
        self.prototype_sets_: tuple[tuple[int, ...], ...] = tuple(
            tuple(sorted(int(index) for index in prototype))
            for prototype in diagnostics["prototypes_original"]
        )
        return self

    def fit_predict(self, X: ArrayLike) -> NDArray[np.int_]:
        """Fit SNI-DPC and return one integer cluster label per sample."""

        return self.fit(X).labels_.copy()


def fit_predict(
    X: ArrayLike,
    *,
    return_diagnostics: bool = False,
) -> NDArray[np.int_] | tuple[NDArray[np.int_], dict[str, Any]]:
    """Cluster a numeric matrix with the frozen SNI-DPC method."""

    estimator = SNIDPC().fit(X)
    labels = estimator.labels_.copy()
    if return_diagnostics:
        return labels, estimator.diagnostics_
    return labels


def _validated_data(X: ArrayLike) -> NDArray[np.float64]:
    data = np.asarray(X, dtype=float)
    if data.ndim != 2:
        raise ValueError("X must be a two-dimensional numeric array.")
    if data.shape[0] < 2:
        raise ValueError("SNI-DPC requires at least two samples.")
    if data.shape[1] < 1:
        raise ValueError("X must contain at least one feature.")
    if not np.isfinite(data).all():
        raise ValueError("X must not contain NaN or infinite values.")
    return np.ascontiguousarray(data, dtype=float)

