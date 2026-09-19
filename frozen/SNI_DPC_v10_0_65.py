"""Frozen SNI-DPC v10_0.65 reference implementation.

This historical reference applies six v10 mechanisms and fixed density weights
to the immutable v9.5 computational core. The reference temporarily overrides
core functions during a call. Use the thread-safe public SNIDPC estimator for
applications. No labels or true cluster count are accepted.
"""

from __future__ import annotations

import math

import numpy as np

from sni_dpc import _core as base


ALGORITHM_VERSION = "SNI-DPC-v10_0.65"


def compute_composite_density(data, NN, RNN, M, S, G):
    """Average the four normalized density views equally."""
    n = len(data)
    rho_nn = np.zeros(n, dtype=float)
    rho_rnn = np.zeros(n, dtype=float)
    rho_snn = np.zeros(n, dtype=float)
    rho_g = np.zeros(n, dtype=float)
    for i in range(n):
        rho_nn[i] = base._distance_density(data, i, NN[i])
        rho_rnn[i] = base._distance_density(data, i, RNN[i])
        rho_snn[i] = base._distance_density(data, i, S[i])
        if M[i]:
            values = [G.get((i, int(j)), 0.0) for j in M[i]]
            rho_g[i] = np.log1p(np.mean(values))
    views = np.vstack(
        [base.minmax_normalize(view) for view in (rho_nn, rho_rnn, rho_snn, rho_g)]
    )
    rho = np.average(views, axis=0, weights=np.asarray([0.35 / 3.0] * 3 + [0.65], dtype=float))
    detail = {
        "rho_nn": rho_nn,
        "rho_rnn": rho_rnn,
        "rho_snn": rho_snn,
        "rho_g": rho_g,
        "rho": rho,
    }
    return rho, detail


def compute_gravity(data, M, sigma, reliability, local_scale):
    """Compute the unit-based reliability-modulated screened interaction."""
    n = len(data)
    G = {}
    G_norm = {}
    mean_G = np.zeros(n, dtype=float)
    for i in range(n):
        values = []
        for j in M[i]:
            j = int(j)
            if i == j:
                continue
            sij = base.get_sigma(sigma, i, j)
            if sij <= base.EPS:
                gij = 0.0
            else:
                dij = base.euclidean_distance(data, i, j)
                scale = math.sqrt(local_scale[i] * local_scale[j]) + base.EPS
                d_scaled = dij / scale
                alpha_ij = base.adaptive_alpha(
                    i, j, reliability, alpha_min=1.0, alpha_max=2.0
                )
                rel_ij = min(reliability[i], reliability[j])
                beta_local = 1.0 + (1.0 - rel_ij)
                gij = (
                    sij**beta_local
                    * math.exp(-(d_scaled**2))
                    / (d_scaled + base.EPS) ** alpha_ij
                )
            G[i, j] = gij
            values.append(gij)
        mean_G[i] = float(np.mean(values)) if values else 0.0
    for (i, j), value in G.items():
        G_norm[i, j] = value / (mean_G[i] + base.EPS)
    return G, G_norm, mean_G


def prototype_force_to_point(
    data,
    gamma,
    rho,
    nb,
    G_norm,
    point_i,
    prototype,
    reliability=None,
    local_scale=None,
    sigma=None,
):
    prototype = list(prototype)
    if not prototype:
        return 0.0
    dists = np.asarray(
        [base.euclidean_distance(data, point_i, int(j)) for j in prototype],
        dtype=float,
    )
    positive = dists[dists > base.EPS]
    proto_scale = float(np.median(positive)) if len(positive) else 1.0
    force = 0.0
    for j, dij in zip(prototype, dists):
        j = int(j)
        gij = G_norm.get((point_i, j))
        if gij is None:
            if reliability is not None and local_scale is not None and sigma is not None:
                sij = base.get_sigma(sigma, point_i, j)
                if sij <= base.EPS:
                    gij = 0.0
                else:
                    scale = (
                        math.sqrt(local_scale[point_i] * local_scale[j]) + base.EPS
                    )
                    d_scaled = dij / scale
                    alpha_ij = base.adaptive_alpha(
                        point_i,
                        j,
                        reliability,
                        alpha_min=1.0,
                        alpha_max=2.0,
                    )
                    rel_ij = min(reliability[point_i], reliability[j])
                    beta_local = 1.0 + (1.0 - rel_ij)
                    gij = (
                        sij**beta_local
                        * math.exp(-(d_scaled**2))
                        / (d_scaled + base.EPS) ** alpha_ij
                    )
            else:
                similarity = base.structure_similarity_from_gamma(
                    gamma, point_i, j
                )
                gij = similarity / (dij * dij + base.EPS)
        proto_wd = math.exp(-(dij * dij) / (proto_scale * proto_scale + base.EPS))
        mass = rho[j] * (1.0 + nb[j])
        if reliability is not None:
            mass *= 1.0 + reliability[j]
        force += mass * gij * proto_wd
    return force / (math.sqrt(len(prototype)) + base.EPS)


def layer1_core_expansion(
    center, labels, current_label, M, S, G_norm, rho, component_ids
):
    center = int(center)
    comp_id = component_ids[center]
    comp_points = np.where(component_ids == comp_id)[0]
    rho_floor = float(np.percentile(rho[comp_points], 25))
    thresholds = np.zeros(len(M), dtype=float)
    for i in range(len(M)):
        values = [
            G_norm.get((i, int(j)), 0.0)
            for j in M[i]
            if G_norm.get((i, int(j)), 0.0) > 0
        ]
        thresholds[i] = float(np.percentile(values, 25)) if values else 0.0
    core = {center}
    frontier = {center}
    while frontier:
        new_nodes = set()
        candidates = set()
        for u in frontier:
            candidates.update(M[u])
        for j in candidates:
            j = int(j)
            if component_ids[j] != comp_id or j in core:
                continue
            if labels[j] != -1 and labels[j] != current_label:
                continue
            strict_hit = any(j in S[u] or u in S[j] for u in core)
            support = max(
                (G_norm.get((j, u), 0.0) for u in core), default=0.0
            )
            if strict_hit or (support >= thresholds[j] and rho[j] >= rho_floor):
                new_nodes.add(j)
        if not new_nodes:
            break
        core.update(new_nodes)
        frontier = new_nodes
    return core


def layer2_prototype_expansion(
    center, core, labels, current_label, data, rho, M, gamma, G_norm, component_ids
):
    center = int(center)
    comp_id = component_ids[center]
    prototype = set(core)
    while True:
        candidates = set()
        for u in prototype:
            candidates.update(M[u])
        second = set()
        for u in candidates:
            second.update(M[int(u)])
        candidates.update(second)
        candidates = {
            int(v)
            for v in candidates
            if component_ids[int(v)] == comp_id
            and v not in prototype
            and (labels[int(v)] == -1 or labels[int(v)] == current_label)
        }
        if not candidates:
            break
        gamma_proto = set()
        for u in prototype:
            gamma_proto.update(gamma[u])
        scores = {}
        for i in candidates:
            values = np.asarray(
                [
                    max(G_norm.get((i, u), 0.0), G_norm.get((u, i), 0.0))
                    for u in prototype
                ],
                dtype=float,
            )
            g_max = float(np.max(values)) if len(values) else 0.0
            g_mean = float(np.mean(values)) if len(values) else 0.0
            similarity = len(gamma[i] & gamma_proto) / (
                len(gamma[i] | gamma_proto) + base.EPS
            )
            ratio = min(1.0, rho[i] / (rho[center] + base.EPS))
            scores[i] = (g_max + g_mean + similarity) / 3.0 * ratio
        threshold = float(np.percentile(list(scores.values()), 25))
        new_nodes = {i for i, score in scores.items() if score >= threshold}
        if not new_nodes:
            break
        old_size = len(prototype)
        prototype.update(new_nodes)
        if len(prototype) == old_size:
            break
    return prototype


def layer3_gravity_boundary_assignment(
    prototype,
    labels,
    current_label,
    prototypes,
    data,
    gamma,
    rho,
    nb,
    M,
    G_norm,
    component_ids,
    reliability=None,
    local_scale=None,
    sigma=None,
):
    prototype = set(prototype)
    if not prototype:
        return set()
    comp_id = component_ids[next(iter(prototype))]
    candidates = set()
    for u in prototype:
        candidates.update(M[int(u)])
    second = set()
    for u in candidates:
        second.update(M[int(u)])
    candidates.update(second)
    candidates = {
        int(i)
        for i in candidates
        if labels[int(i)] == -1
        and component_ids[int(i)] == comp_id
        and int(i) not in prototype
    }
    if not candidates:
        return set()
    same = [
        item
        for item in prototypes
        if item and component_ids[next(iter(item))] == comp_id
    ]
    all_prototypes = same + [prototype]
    current_index = len(all_prototypes) - 1
    current = []
    for i in candidates:
        force = prototype_force_to_point(
            data,
            gamma,
            rho,
            nb,
            G_norm,
            i,
            prototype,
            reliability=reliability,
            local_scale=local_scale,
            sigma=sigma,
        )
        if force <= base.EPS:
            distance = min(
                base.euclidean_distance(data, i, int(p)) for p in prototype
            )
            force = 1.0 / (distance + base.EPS)
        current.append(force)
    positive = np.asarray([value for value in current if value > base.EPS])
    if not len(positive):
        return set()
    threshold = float(np.percentile(positive, 25))
    accepted = set()
    for i in candidates:
        forces = []
        for item in all_prototypes:
            force = prototype_force_to_point(
                data,
                gamma,
                rho,
                nb,
                G_norm,
                i,
                item,
                reliability=reliability,
                local_scale=local_scale,
                sigma=sigma,
            )
            if force <= base.EPS:
                distance = min(
                    base.euclidean_distance(data, i, int(p)) for p in item
                )
                force = 1.0 / (distance + base.EPS)
            forces.append(force)
        if int(np.argmax(forces)) == current_index and forces[current_index] >= threshold:
            accepted.add(i)
    return accepted


def mng_dpc(data):
    """Run frozen SNI-DPC v10_0.65 without editing the historical core."""
    overrides = {
        "compute_composite_density": compute_composite_density,
        "compute_gravity": compute_gravity,
        "prototype_force_to_point": prototype_force_to_point,
        "layer1_core_expansion": layer1_core_expansion,
        "layer2_prototype_expansion": layer2_prototype_expansion,
        "layer3_gravity_boundary_assignment": layer3_gravity_boundary_assignment,
    }
    original = {name: getattr(base, name) for name in overrides}
    try:
        for name, function in overrides.items():
            setattr(base, name, function)
        return base.mng_dpc(data)
    finally:
        for name, function in original.items():
            setattr(base, name, function)
