import numpy as np
from scipy.spatial import KDTree


def build_kdtree(events):
    """
    Build a KDTree for fast nearest-neighbour lookup.

    Rapidity (y) is scaled so that the distance metric matches:
    d = sqrt((pt_i - pt_j)^2 + 100*(y_i - y_j)^2)

    Scaling y by 10 makes Euclidean distance equivalent.
    """
    points = np.column_stack((events["pt"], 10 * events["y"]))
    tree = KDTree(points)
    return tree


def redistribute_weights(cell_indices, events):
    """
    Apply the weight redistribution formula:

    w'_i = |w_i| / sum_j(|w_j|) * sum_j(w_j)

    This guarantees:
    - all weights become positive
    - total weight is conserved
    """

    weights = events.loc[cell_indices, "w"].values

    S = np.sum(weights)
    S_abs = np.sum(np.abs(weights))

    new_weights = np.abs(weights) / S_abs * S

    events.loc[cell_indices, "w"] = new_weights


def grow_cell(seed_index, events, tree, used_mask):
    """
    Grow a cell around a negative-weight seed event.

    Keep adding nearest neighbours until the total weight
    inside the cell becomes >= 0.

    Uses incremental k-doubling to avoid querying all neighbours
    upfront in the common case.
    """

    seed_point = np.array([events.loc[seed_index, "pt"],
                           10 * events.loc[seed_index, "y"]])

    n_total = len(events)
    k = min(10, n_total)
    fetched_indices = []

    while True:
        _, indices = tree.query(seed_point, k=k)
        fetched_indices = indices if np.ndim(indices) == 1 else indices.ravel()

        cell = []
        weight_sum = 0.0
        found = False

        for idx in fetched_indices:
            if used_mask[idx]:
                continue
            cell.append(idx)
            weight_sum += events.loc[idx, "w"]
            if weight_sum >= 0:
                found = True
                break

        if found:
            return cell

        if k >= n_total:
            # All events exhausted without reaching non-negative sum
            raise ValueError(
                f"Cell starting at seed index {seed_index} could not reach "
                f"non-negative weight sum (total weight = {weight_sum:.4f}). "
                "The global weight sum may be negative."
            )

        k = min(k * 2, n_total)


def cell_resampling(events):
    """
    Main algorithm.

    Input:
        DataFrame with columns: pt, y, weight

    Output:
        DataFrame with resampled weights (all >= 0)
    """

    events = events.copy().reset_index(drop=True)
    events = events.rename(columns={"weight": "w"})

    tree = build_kdtree(events)

    used_mask = np.zeros(len(events), dtype=bool)

    for i in range(len(events)):

        if used_mask[i]:
            continue

        if events.loc[i, "w"] >= 0:
            continue

        # grow a cell starting from negative seed
        cell_indices = grow_cell(i, events, tree, used_mask)

        # redistribute weights
        redistribute_weights(cell_indices, events)

        # mark events as used
        for idx in cell_indices:
            used_mask[idx] = True

    events = events.rename(columns={"w": "weight"})
    return events