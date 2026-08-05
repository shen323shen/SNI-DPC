# -*- coding: utf-8 -*-
"""
SNI-DPC v9.5: Screened Natural-Neighbor Interaction Density Peaks Clustering.

Publication copy of the frozen experimental implementation. The explicit
UTF-8 declaration prevents legacy Windows tools from displaying the existing
Chinese implementation notes with the system code page. Computational logic
is unchanged from experiments/src/MNG_DPC_v9_5.py.

设计目标：
1. 不输入簇个数；
2. 不设置人为开关；
3. 按模块组织代码；
4. 数据集路径放在文件最后；
5. 保留 DBCP / DSNGCAP 的主线：NaNN -> 密度低置信点识别 -> clean 图重建 -> 聚类 -> 低密度点重分配；
6. 将 getSNG 的普通连通分量扩展，替换为：中心探测 + 第一层簇核 + 第二层簇原型 + 第三层引力边界分配。

依赖：
    numpy scipy scikit-learn matplotlib pandas

说明：
    这是第一版“逻辑闭环代码”，优先保证算法流程完整、变量清晰、便于后续逐模块调参/消融。
"""

import os
import math
import time
from collections import deque, Counter

import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA


EPS = 1e-12
ALGORITHM_VERSION = "9.5"


from collections import defaultdict, Counter


def component_shape_info(data, comp):
    idx = np.asarray(list(comp), dtype=int)

    if len(idx) <= 3:
        return {
            "thin_ratio": 1.0,
            "axis": np.ones(data.shape[1], dtype=float),
            "center": np.mean(data[idx], axis=0),
        }

    X = data[idx] - np.mean(data[idx], axis=0)

    try:
        _, s, vt = np.linalg.svd(X, full_matrices=False)

        if len(s) < 2 or s[0] <= EPS:
            thin_ratio = 0.0
        else:
            thin_ratio = float(s[1] / (s[0] + EPS))

        axis = vt[0]
        axis = axis / (np.linalg.norm(axis) + EPS)

    except Exception:
        thin_ratio = 1.0
        axis = np.ones(data.shape[1], dtype=float)
        axis = axis / (np.linalg.norm(axis) + EPS)

    return {
        "thin_ratio": thin_ratio,
        "axis": axis,
        "center": np.mean(data[idx], axis=0),
    }

def build_raw_to_cut_component_map(component_ids_raw, components_cut):
    raw_to_cut = {}

    for cut_cid, comp in enumerate(components_cut):
        raw_ids = [int(component_ids_raw[int(i)]) for i in comp]
        raw_cid = Counter(raw_ids).most_common(1)[0][0]

        if raw_cid not in raw_to_cut:
            raw_to_cut[raw_cid] = []

        raw_to_cut[raw_cid].append(cut_cid)

    return raw_to_cut


def min_inter_component_distance(data, comp_a, comp_b):
    idx_a = np.asarray(list(comp_a), dtype=int)
    idx_b = np.asarray(list(comp_b), dtype=int)

    if len(idx_a) == 0 or len(idx_b) == 0:
        return np.inf

    tree = KDTree(data[idx_b])
    dist, _ = tree.query(data[idx_a], k=1)

    return float(np.min(dist))


def axis_similarity(axis_a, axis_b):
    axis_a = np.asarray(axis_a, dtype=float)
    axis_b = np.asarray(axis_b, dtype=float)

    return float(
        abs(np.dot(axis_a, axis_b))
        / ((np.linalg.norm(axis_a) + EPS) * (np.linalg.norm(axis_b) + EPS))
    )


def compute_raw_component_center_caps_v94(
    data,
    components_raw,
    components_cut,
    component_ids_raw,
    ref_size,
    local_scale=None
):
    """
    v9.4 中心预算：

    1. raw component 仍然决定中心预算上限；
    2. 但不再只按 size / ref_size 给中心；
    3. 如果 raw component 是细长线型结构，则强制 cap = 1；
    4. 如果 raw component 只有一个显著 cut component，则 cap = 1；
    5. 很小且贴近同方向大线型 component 的 raw component，不单独给中心。
    """
    raw_to_cut = build_raw_to_cut_component_map(
        component_ids_raw,
        components_cut
    )

    raw_shapes = {}
    raw_center_caps = {}
    raw_parent = {}

    # 先初始化 parent
    for cid in range(len(components_raw)):
        raw_parent[cid] = cid
        raw_shapes[cid] = component_shape_info(data, components_raw[cid])

    # 自适应的“显著 cut component”阈值
    significant_cut_min = max(25.0, 0.50 * float(ref_size))

    for raw_cid, comp in enumerate(components_raw):
        size = len(comp)
        shape = raw_shapes[raw_cid]

        size_cap = int(round(size / (ref_size + EPS)))
        size_cap = max(1, size_cap)

        cut_ids = raw_to_cut.get(raw_cid, [])
        cut_sizes = [len(components_cut[cid]) for cid in cut_ids]

        significant_cut_count = sum(
            1 for s in cut_sizes
            if s >= significant_cut_min
        )

        significant_cut_count = max(1, significant_cut_count)

        # 先用 size_cap 和显著 cut 数共同约束
        cap = min(size_cap, significant_cut_count)

        # 线型保护：
        # 细长 component 即使样本很多，也更可能是一个长条簇，而不是多个簇。
        if shape["thin_ratio"] <= 0.08 and size >= 1.20 * ref_size:
            cap = 1

        raw_center_caps[raw_cid] = cap

    # 小 raw component 吸收：
    # 如果一个很小的 raw component 靠近某个大线型 component，
    # 且方向基本一致，则认为它是该线型簇的碎片，不单独给中心。
    main_raw_ids = [
        cid for cid, comp in enumerate(components_raw)
        if len(comp) >= 0.75 * ref_size
    ]

    for raw_cid, comp in enumerate(components_raw):
        size = len(comp)

        if size >= 0.45 * ref_size:
            continue

        if len(main_raw_ids) == 0:
            continue

        best_parent = None
        best_dist = np.inf

        for pid in main_raw_ids:
            if pid == raw_cid:
                continue

            d = min_inter_component_distance(
                data,
                components_raw[raw_cid],
                components_raw[pid]
            )

            if d < best_dist:
                best_dist = d
                best_parent = pid

        if best_parent is None:
            continue

        shape_small = raw_shapes[raw_cid]
        shape_big = raw_shapes[best_parent]

        same_axis = axis_similarity(
            shape_small["axis"],
            shape_big["axis"]
        ) >= 0.90

        both_line_like = (
            shape_small["thin_ratio"] <= 0.12
            and shape_big["thin_ratio"] <= 0.12
        )

        if local_scale is not None:
            idx_small = np.asarray(list(comp), dtype=int)
            local_gap_scale = float(np.median(local_scale[idx_small]))
        else:
            local_gap_scale = 1.0

        near_parent = best_dist <= 3.5 * (local_gap_scale + EPS)

        if both_line_like and same_axis and near_parent:
            raw_center_caps[raw_cid] = 0
            raw_parent[raw_cid] = best_parent

    cap_info = {
        "raw_shapes": raw_shapes,
        "raw_to_cut": raw_to_cut,
        "raw_parent": raw_parent,
        "significant_cut_min": significant_cut_min,
    }

    return raw_center_caps, cap_info




def select_centers_by_raw_cut_budget(
    base_score,
    components_cut,
    component_ids_raw,
    components_raw,
    raw_center_caps
):
    """
    v9.3 中心选择策略：

    1. raw component 决定最多几个中心；
    2. M_cc component 决定中心候选区域；
    3. 每个 raw component 内，优先从最大的 cut components 中选中心；
    4. 不让 tiny cut component 都变成中心。
    """
    raw_to_cut = build_raw_to_cut_component_map(
        component_ids_raw,
        components_cut
    )

    center_records = []

    for raw_cid in range(len(components_raw)):
        cut_ids = raw_to_cut.get(raw_cid, [])

        if len(cut_ids) == 0:
            continue

        # 优先选择大 cut component，避免小碎片产生中心
        cut_ids = sorted(
            cut_ids,
            key=lambda cid: len(components_cut[cid]),
            reverse=True
        )

        quota = int(raw_center_caps.get(raw_cid, 1))
        if quota == 0:
            continue
        quota = min(quota, len(cut_ids))

        for cut_cid in cut_ids[:quota]:
            idx = np.asarray(list(components_cut[cut_cid]), dtype=int)

            if len(idx) == 0:
                continue

            local_best = int(idx[np.argmax(base_score[idx])])

            center_records.append({
                "center": local_best,
                "raw_cid": raw_cid,
                "cut_cid": cut_cid,
                "cut_size": len(idx),
                "score": float(base_score[local_best]),
            })

    # 高分中心先扩展，但所有中心会提前保护，所以顺序影响会小很多
    center_records = sorted(
        center_records,
        key=lambda r: r["score"],
        reverse=True
    )

    centers = [r["center"] for r in center_records]

    return centers, center_records


# ============================================================
# 1. 数据读取
# ============================================================

def load_data(data_path: str):
    """读取 csv / txt / mat 数据。默认前两列或 mat 中 data 字段为特征。"""
    suffix = data_path.split(".")[-1].lower()

    if suffix in ["csv", "txt", "data"]:
        try:
            data = np.asarray(np.loadtxt(data_path, delimiter=",", comments="x", usecols=[0, 1]), dtype=float)
        except Exception:
            df = pd.read_csv(data_path, header=None)
            data = df.iloc[:, :2].values.astype(float)
    elif suffix == "mat":
        data_dict = sio.loadmat(data_path)
        if "data" in data_dict:
            data = data_dict["data"]
        elif "X" in data_dict:
            data = data_dict["X"]
        else:
            raise KeyError("mat 文件中没有找到 data 或 X 字段，请手动修改 load_data。")
        data = np.asarray(data, dtype=float)
    else:
        raise ValueError(f"暂不支持该文件格式: {suffix}")

    if data.ndim != 2:
        raise ValueError("data 必须是二维数组。")
    return data


# ============================================================
# 2. 基础工具函数
# ============================================================

def safe_set(values):
    return set(int(v) for v in values if v is not None)


def minmax_normalize(arr):
    arr = np.asarray(arr, dtype=float)
    mn, mx = np.min(arr), np.max(arr)
    if abs(mx - mn) < EPS:
        return np.zeros_like(arr, dtype=float)
    return (arr - mn) / (mx - mn + EPS)


def euclidean_distance(data, i, j):
    return float(np.linalg.norm(data[i] - data[j]))


def structure_similarity_from_gamma(gamma, i, j, mode="scan"):
    """
    结构相似性。
    默认使用 SCAN 型结构相似性：|Γ_i ∩ Γ_j| / sqrt(|Γ_i||Γ_j|)。
    mode 参数内部固定默认，不在主流程中作为开关使用。
    """
    gi, gj = gamma[i], gamma[j]
    if len(gi) == 0 or len(gj) == 0:
        return 0.0

    inter = len(gi & gj)
    if mode == "jaccard":
        union = len(gi | gj)
        return inter / (union + EPS)
    return inter / (math.sqrt(len(gi) * len(gj)) + EPS)


def auto_upper_threshold(values):
    """
    根据排序曲线最大间隔自动取高分候选阈值。
    用于中心分数，不输入簇个数。
    """
    values = np.asarray(values, dtype=float)
    if len(values) <= 2:
        return np.max(values)

    sorted_values = np.sort(values)
    diffs = np.diff(sorted_values)
    if np.max(diffs) < EPS:
        return np.max(values) + 1.0

    # 避免极低端的小波动影响，只在后 70% 范围找最大跳跃，更偏向挑出高分中心候选。
    start = max(0, int(len(diffs) * 0.30))
    pos = start + int(np.argmax(diffs[start:]))
    threshold = sorted_values[pos + 1]
    return threshold


# ============================================================
# 3. 自然邻居搜索 NaNN
# ============================================================

def _query_knn_excluding_self(data, tree, k_neighbors):
    """查询每个点的前 k_neighbors 个近邻，不含自身。"""
    n = len(data)
    k_query = min(n, k_neighbors + 1)
    distances, indices = tree.query(data, k=k_query)

    knn_indices = []
    knn_distances = []
    for i in range(n):
        row_idx = []
        row_dist = []
        for d, j in zip(distances[i], indices[i]):
            j = int(j)
            if j == i:
                continue
            row_idx.append(j)
            row_dist.append(float(d))
            if len(row_idx) >= k_neighbors:
                break
        knn_indices.append(row_idx)
        knn_distances.append(row_dist)
    return knn_indices, knn_distances


def natural_neighbor_search(data):
    """
    最近自然邻居搜索。

    返回：
        r: 自然稳定轮次
        NN: 每个点主动选择的自然邻居集合 list[set]
        RNN: 每个点的逆邻居集合 list[set]
        nb: 每个点逆邻居数量 np.ndarray
        tree: KDTree
    """
    data = np.asarray(data, dtype=float)
    n = len(data)
    tree = KDTree(data)

    if n <= 1:
        return 0, [set() for _ in range(n)], [set() for _ in range(n)], np.zeros(n, dtype=int), tree

    # 初始查询上限。若自然稳定轮次超过该值，会自动扩展。
    kmax = min(n - 1, max(8, int(math.sqrt(n) * 4)))
    knn_indices, _ = _query_knn_excluding_self(data, tree, kmax)

    NN = [set() for _ in range(n)]
    RNN = [set() for _ in range(n)]
    zero_history = []

    r = 0
    while True:
        r += 1
        if r > kmax:
            if kmax >= n - 1:
                break
            kmax = min(n - 1, kmax * 2)
            knn_indices, _ = _query_knn_excluding_self(data, tree, kmax)

        for i in range(n):
            if r - 1 < len(knn_indices[i]):
                j = int(knn_indices[i][r - 1])
                NN[i].add(j)
                RNN[j].add(i)

        nb = np.asarray([len(RNN[i]) for i in range(n)], dtype=int)
        zero_count = int(np.sum(nb == 0))

        if zero_count == 0:
            break

        zero_history.append(zero_count)
        # 与参考代码一致：连续三轮 0-nb 点数不再变化则停止。
        if len(zero_history) >= 3 and zero_history[-1] == zero_history[-2] == zero_history[-3]:
            break

        if r >= n - 1:
            break

    nb = np.asarray([len(RNN[i]) for i in range(n)], dtype=int)
    return r, NN, RNN, nb, tree


# ============================================================
# 4. 多邻域图 MNG 构建
# ============================================================

def build_mng(NN, RNN):
    """
    构建多邻域图。

    M_i = NN_i ∪ RNN_i
    S_i = NN_i ∩ RNN_i，即严格自然邻居集合。
    Γ_i = {i} ∪ M_i，用于结构相似性。
    """
    n = len(NN)
    M = []
    S = []
    gamma = []

    for i in range(n):
        nn_i = safe_set(NN[i])
        rnn_i = safe_set(RNN[i])
        m_i = nn_i | rnn_i
        s_i = nn_i & rnn_i
        M.append(m_i)
        S.append(s_i)
        gamma.append(set([i]) | m_i)

    return M, S, gamma


# ============================================================
# 5. 结构相似性与引力
# ============================================================
def compute_reliability_from_nb(nb, eps=EPS):
    """
    根据逆邻居数量 nb 构造点可靠性权重。

    nb_i 越大，说明点 i 被更多点当作邻居，
    越像局部核心点，可靠性越高。

    reliability_i ∈ (0, 1]
    """
    nb = np.asarray(nb, dtype=float)
    nb_log = np.log1p(nb)

    max_v = np.max(nb_log)
    if max_v <= eps:
        return np.ones_like(nb_log)

    reliability = nb_log / (max_v + eps)

    # 防止 0 权重
    reliability = np.maximum(reliability, eps)

    return reliability

def compute_local_scale_from_gamma(data, gamma, eps=EPS):
    """
    计算每个点的局部距离尺度 s_i。

    s_i = mean_{u in Gamma_i, u != i} d(i,u)

    作用：
    - 密集区域 s_i 小；
    - 稀疏区域 s_i 大；
    - 后续用 d_ij / sqrt(s_i*s_j) 做局部尺度归一化。
    """
    # 局部自适应距离度量的关键步骤
    n = len(data)
    local_scale = np.zeros(n, dtype=float)
    all_scales = []

    for i in range(n):
        dists = []

        for u in gamma[i]:
            u = int(u)
            if u == i:
                continue

            d = euclidean_distance(data, i, u)
            if d > eps:
                dists.append(d)

        if len(dists) > 0:
            local_scale[i] = float(np.mean(dists))
            all_scales.append(local_scale[i])
        else:
            local_scale[i] = 0.0

    if len(all_scales) > 0:
        fallback = float(np.median(all_scales))
    else:
        fallback = 1.0

    fallback = max(fallback, eps)

    local_scale[local_scale <= eps] = fallback

    return local_scale

def weighted_common_neighbor_similarity(gamma_i, gamma_j, reliability, eps=EPS):
    """
    可靠共同邻居结构相似度。

    原始：
        sigma_ij = |Gamma_i ∩ Gamma_j| / sqrt(|Gamma_i| |Gamma_j|)

    改进：
        RCN_ij = sum_{u in common} w_u
                 / sqrt(sum_{u in Gamma_i} w_u * sum_{u in Gamma_j} w_u)
    """
    common = gamma_i & gamma_j

    if len(common) == 0:
        return 0.0

    common_weight = sum(reliability[int(u)] for u in common)

    weight_i = sum(reliability[int(u)] for u in gamma_i)
    weight_j = sum(reliability[int(u)] for u in gamma_j)

    return float(common_weight / (math.sqrt(weight_i * weight_j) + eps))

def adaptive_alpha(i, j, reliability, alpha_min=1.3, alpha_max=2.5):
    """
    自适应距离衰减指数。

    reliability 高：更像核心点，距离惩罚稍弱；
    reliability 低：更像边界点或噪声点，距离惩罚更强。
    """
    ri = reliability[int(i)]
    rj = reliability[int(j)]

    alpha_i = alpha_min + (1.0 - ri) * (alpha_max - alpha_min)
    alpha_j = alpha_min + (1.0 - rj) * (alpha_max - alpha_min)

    return float((alpha_i + alpha_j) / 2.0)






def compute_structure_similarity(M, gamma, reliability):
    """
    只对 MNG 中出现的边计算结构相似性。

    这里不再使用普通共同邻居数量，
    而是使用“可靠共同邻居加权结构相似度”。
    """
    sigma = {}
    n = len(M)

    for i in range(n):
        for j in M[i]:
            j = int(j)

            if i == j:
                continue

            key = (min(i, j), max(i, j))

            if key not in sigma:
                sigma[key] = weighted_common_neighbor_similarity(
                    gamma[i],
                    gamma[j],
                    reliability
                )

    return sigma

def get_sigma(sigma, i, j):
    if i == j:
        return 1.0
    return float(sigma.get((min(i, j), max(i, j)), 0.0))


def compute_gravity(
    data,
    M,
    sigma,
    reliability,
    local_scale,
    beta=1.5,
    alpha_min=1.3,
    alpha_max=2.5
):
    """Compute the screened interaction kernel on directed MNG edges.

    K_ij = sigma_ij**beta_ij * exp(-z_ij**2)
           / (z_ij + EPS)**alpha_ij,
    where z_ij is the distance scaled by the geometric mean of the two
    endpoint-local scales.  No separate excess-distance barrier is used.
    """
    n = len(data)
    G = {}
    G_norm = {}
    mean_G = np.zeros(n, dtype=float)

    for i in range(n):
        row_values = []

        for j in M[i]:
            j = int(j)

            if i == j:
                continue

            sij = get_sigma(sigma, i, j)

            if sij <= EPS:
                gij = 0.0
                G[(i, j)] = gij
                row_values.append(gij)
                continue

            dij = euclidean_distance(data, i, j)

            pair_scale = math.sqrt(local_scale[i] * local_scale[j]) + EPS
            d_scaled = dij / pair_scale

            alpha_ij = adaptive_alpha(
                i,
                j,
                reliability,
                alpha_min=alpha_min,
                alpha_max=alpha_max
            )

            # Smooth Gaussian screening suppresses long edges without a threshold.
            wd = math.exp(-(d_scaled ** 2))

            # Endpoint reliability adapts beta so weak boundary edges receive
            # a stronger structural-support penalty.
            rel_ij = min(reliability[i], reliability[j])
            beta_local = beta + 0.3 + (1.0 - rel_ij) * 0.7

            # Adaptive beta down-weights edges with weak structural support.
            gij = (sij ** beta_local) * wd / ((d_scaled + EPS) ** alpha_ij)

            G[(i, j)] = gij
            row_values.append(gij)

        if len(row_values) > 0:
            mean_G[i] = float(np.mean(row_values))
        else:
            mean_G[i] = 0.0

    # 出边均值归一化
    for (i, j), gij in G.items():
        G_norm[(i, j)] = gij / (mean_G[i] + EPS)

    return G, G_norm, mean_G


# ============================================================
# 6. 复合密度
# ============================================================

def _distance_density(data, i, neighbors):
    neighbors = list(neighbors)
    if len(neighbors) == 0:
        return 0.0
    dist_sum = 0.0
    for j in neighbors:
        dist_sum += euclidean_distance(data, i, j)
    return len(neighbors) / (dist_sum + EPS)


def compute_composite_density(data, NN, RNN, M, S, G):
    """Combine three graph densities with a stable interaction density.

    The interaction term is log(1 + mean_j K_ij), avoiding the algebraic
    saturation caused by averaging row-normalized interactions.
    """
    n = len(data)
    rho_nn = np.zeros(n, dtype=float)
    rho_rnn = np.zeros(n, dtype=float)
    rho_snn = np.zeros(n, dtype=float)
    rho_g = np.zeros(n, dtype=float)

    for i in range(n):
        rho_nn[i] = _distance_density(data, i, NN[i])
        rho_rnn[i] = _distance_density(data, i, RNN[i])
        rho_snn[i] = _distance_density(data, i, S[i])
        if len(M[i]) > 0:
            raw_mean = np.mean([G.get((i, j), 0.0) for j in M[i]])
            rho_g[i] = math.log1p(raw_mean)
        else:
            rho_g[i] = 0.0

    rho = (
        0.20 * minmax_normalize(rho_nn)
        + 0.20 * minmax_normalize(rho_rnn)
        + 0.20 * minmax_normalize(rho_snn)
        + 0.40 * minmax_normalize(rho_g)
    )

    detail = {
        "rho_nn": rho_nn,
        "rho_rnn": rho_rnn,
        "rho_snn": rho_snn,
        "rho_g": rho_g,
        "rho": rho,
    }
    return rho, detail



# ============================================================
# 7. 低密度 / 低置信点自动识别
# ============================================================
#DSNGCAP 密度计算
def compute_dsngcap_density(data, tree, N):
    """
    DSNGCAP 密度：
        density_i = N / sum_{j in N_i} d_ij

    N 建议取：
        N = (r + max(nb)) // 2
    """
    data = np.asarray(data, dtype=float)
    n = len(data)

    N = int(max(1, min(N, n - 1)))
    density = np.zeros(n, dtype=float)

    for i in range(n):
        dist, idx = tree.query([data[i]], k=N + 1)
        dist_list = dist[0][1:]   # 去掉自己

        dist_sum = np.sum(dist_list)
        if dist_sum <= EPS:
            density[i] = 1e8
        else:
            density[i] = N / (dist_sum + EPS)

    return density
#一阶差分
def compute_diff1_by_density(density):
    """
    对密度从小到大排序后，计算一阶差分。
    """
    sorted_density = np.sort(density)
    diff1 = []

    for i in range(len(sorted_density) - 1):
        diff1.append(sorted_density[i + 1] - sorted_density[i])

    return np.asarray(diff1, dtype=float)
#分段方差
def dsngcap_var_compute(diff1, scale):
    """
    将一阶差分 diff1 分成 scale 段，计算每一段的方差。
    """
    diff1 = list(diff1)

    if len(diff1) == 0:
        return np.asarray([0.0])

    diff1.append(diff1[-1])

    step = max(1, len(diff1) // scale)
    var = []

    for i in range(scale):
        start = i * step
        end = min((i + 1) * step, len(diff1))

        if start >= len(diff1):
            var.append(0.0)
        else:
            temp_diff1 = diff1[start:end]
            var.append(float(np.var(temp_diff1)))

    return np.asarray(var, dtype=float)
#DSNGCAP 的 alpha 寻找逻辑
def dsngcap_get_alpha(var, scale):
    """
    DSNGCAP 的密度跳变点检测。
    返回低密度点比例 alpha，单位：百分比。
    """
    if len(var) == 0:
        return 0

    standard = scale // 2
    var = list(var[:standard])

    if len(var) == 0:
        return 0

    if max(var) <= EPS:
        return 0

    # 1. 找前半段方差波峰
    max_index = int(np.argmax(var))

    # 2. 从波峰后寻找稳定区
    var_range = max(var) - var[-1]
    k = max(1, standard - max_index)

    over = []

    for i in range(1, k + 1):
        step_value = (var_range / k) * i
        p = 0

        for j in range(standard - 1, max_index - 1, -1):
            if step_value > var[j]:
                p += 1
            else:
                break

        over.append(p)

    if len(over) == 0:
        return 0

    enhance = []
    enhance.append(over[0])

    for i in range(1, len(over)):
        if over[i] == 0 or over[i - 1] == 0:
            enhance.append(1.0)
        else:
            enhance.append(over[i] / (over[i - 1] + EPS))

    a = int(np.argmax(enhance))
    alpha = standard - over[a]

    if max_index != 0:
        alpha -= 1

    alpha = max(0, alpha)

    if scale == 100:
        return int(alpha)
    else:
        return int(alpha * 10)

#最终 DSNGCAP 去噪函数
def detect_low_density_points_by_dsngcap_density(data, tree, r, nb):
    """
    直接采用 DSNGCAP 的密度跳变点去噪。
    不使用 Q，不使用原始置信度，不使用 rho * nb * mean_G。
    """
    n = len(data)

    if n <= 5:
        return np.asarray([], dtype=int), np.arange(n), 0, None

    # DSNGCAP 中常用的 N
    N = int((r + int(np.max(nb))) // 2)
    N = max(1, min(N, n - 1))

    # 1. 计算 DSNGCAP 密度
    density = compute_dsngcap_density(data, tree, N)

    # 2. 密度从小到大排序后计算一阶差分
    diff1 = compute_diff1_by_density(density)

    # 3. 根据数据规模确定尺度
    if n >= 500:
        scale = 100
    else:
        scale = 10

    # 4. 分段方差
    var = dsngcap_var_compute(diff1, scale)

    # 5. 自动寻找低密度比例 alpha
    alpha = dsngcap_get_alpha(var, scale)

    # 6. 根据 alpha 划分低密度点
    sorted_idx = np.argsort(density)
    low_len = int(n * alpha / 100)

    # 防止过度删空
    low_len = max(0, min(low_len, n - 1))

    low_indices = np.asarray(sorted_idx[:low_len], dtype=int)
    clean_indices = np.asarray(sorted_idx[low_len:], dtype=int)

    denoise_info = {
        "N": N,
        "density": density,
        "diff1": diff1,
        "var": var,
        "scale": scale,
        "alpha": alpha,
    }

    return low_indices, clean_indices, alpha, denoise_info


def detect_low_density_points_by_nb_adaptive(nb, r, eps=EPS):
    """
    Gravity 改进版：基于逆邻居数 nb 自适应去噪。

    核心思想（来自 Gravity 的 NoiseSeparate）：
        nb_i（逆邻居数）越低 → 越可能是噪声/离群点。

    改进：
        k 不再人工设定，而是根据 nb 的正值分布自动确定：
        k_auto = max(1, 15%分位数 of nb_positive)
        同时限制最多去除 30% 的点，防止过度去噪。
    """
    nb = np.asarray(nb, dtype=float)
    n = len(nb)

    if n <= 5:
        return np.asarray([], dtype=int), np.arange(n), 0, {"k_auto": 0}

    # 只在 nb>0 的点中取 15% 分位作为自适应阈值
    nb_positive = nb[nb > eps]
    if len(nb_positive) > 0:
        k_auto = max(1, int(np.percentile(nb_positive, 15)))
    else:
        k_auto = 0

    # 标记低逆邻居点
    noise_mask = nb <= k_auto
    low_indices = np.where(noise_mask)[0].astype(int)

    # 防止过度去噪：最多去掉 30%
    max_remove = int(n * 0.30)
    if len(low_indices) > max_remove:
        sorted_idx = np.argsort(nb)
        low_indices = sorted_idx[:max_remove].astype(int)

    clean_indices = np.setdiff1d(np.arange(n, dtype=int), low_indices).astype(int)
    alpha = int(round(len(low_indices) / n * 100))

    denoise_info = {
        "k_auto": k_auto,
        "method": "nb_adaptive",
    }

    return low_indices, clean_indices, alpha, denoise_info




def compute_low_confidence_score(rho, nb, mean_G):
    """
    Q_i = rho_i * (1 + nb_i) * (1 + mean_G_i)
    Q 越小越可能是低密度、低置信点。
    """
    rho_n = minmax_normalize(rho)
    nb_n = minmax_normalize(nb)
    g_n = minmax_normalize(mean_G)
    return rho_n * (1.0 + nb_n) * (1.0 + g_n)


def _var_compute(diff1, scale):
    diff1 = list(diff1)
    if len(diff1) == 0:
        return [0.0]
    diff1.append(diff1[-1])
    step = max(1, len(diff1) // scale)
    var = []
    for i in range(scale):
        start, end = i * step, min((i + 1) * step, len(diff1))
        if start >= len(diff1):
            var.append(0.0)
        else:
            var.append(float(np.var(diff1[start:end])))
    return var


def _get_alpha_from_var(var, scale):
    """
    参考 DSNGCAP 中 getAlpha1 的自动分界思想。
    返回低置信点比例 alpha，单位是百分比。
    """
    if len(var) == 0:
        return 0

    standard = scale // 2
    var = list(var[:standard])
    if len(var) == 0 or max(var) < EPS:
        return 0

    max_index = int(np.argmax(var))
    var_range = max(var) - var[-1]
    k = max(1, standard - max_index)

    over = []
    for i in range(1, k + 1):
        step_value = (var_range / k) * i
        p = 0
        for j in range(standard - 1, max_index - 1, -1):
            if step_value > var[j]:
                p += 1
            else:
                break
        over.append(p)

    if len(over) == 0:
        return 0

    enhance = [over[0]]
    for i in range(1, len(over)):
        if over[i] == 0 or over[i - 1] == 0:
            enhance.append(1.0)
        else:
            enhance.append(over[i] / (over[i - 1] + EPS))

    a = int(np.argmax(enhance))
    alpha = standard - over[a]
    if max_index != 0:
        alpha -= 1

    alpha = max(0, alpha)
    if scale == 100:
        return int(alpha)
    return int(alpha * 10)


def detect_low_density_points(Q):
    """
    按 Q 从小到大排序，用差分方差自动确定低置信点比例。

    返回：
        low_indices: 低置信点原索引
        clean_indices: clean 点原索引
        alpha: 低置信点比例
    """
    n = len(Q)
    if n <= 5:
        return np.asarray([], dtype=int), np.arange(n), 0

    scale = 100 if n >= 500 else 10
    sorted_idx = np.argsort(Q)
    sorted_Q = Q[sorted_idx]
    diff1 = np.diff(sorted_Q)
    var = _var_compute(diff1, scale)
    alpha = _get_alpha_from_var(var, scale)

    low_len = int(n * alpha / 100)
    low_len = max(0, min(low_len, n - 1))

    low_indices = np.asarray(sorted_idx[:low_len], dtype=int)
    clean_indices = np.asarray(sorted_idx[low_len:], dtype=int)
    return low_indices, clean_indices, alpha


# ============================================================
# 8. 依赖向量、相对距离、中心得分
# ============================================================
#新增弱连通分量函数
def weak_connected_components(M):
    """
    将 MNG 有向图转成弱连通图，计算每个点所属连通分量。

    返回：
        component_ids: 每个点所属连通分量编号
        components: 每个连通分量包含的点集合
    """
    n = len(M)

    undirected = [set() for _ in range(n)]

    for i in range(n):
        for j in M[i]:
            j = int(j)
            if i == j:
                continue
            undirected[i].add(j)
            undirected[j].add(i)

    component_ids = np.full(n, -1, dtype=int)
    components = []
    comp_id = 0

    for start in range(n):
        if component_ids[start] != -1:
            continue

        q = deque([start])
        component_ids[start] = comp_id
        comp = set([start])

        while q:
            i = q.popleft()

            for j in undirected[i]:
                if component_ids[j] == -1:
                    component_ids[j] = comp_id
                    comp.add(j)
                    q.append(j)

        components.append(comp)
        comp_id += 1

    return component_ids, components

def robust_reference_component_size(components_cut):
    """
    从剪桥后的 M_cc components 中估计典型真实簇规模。

    思路：
    - M_cc 可能产生很多小碎片；
    - 小碎片不能参与典型簇大小估计；
    - 因此取 component size 的上半部分，再取中位数。
    """
    if len(components_cut) == 0:
        return 1.0

    sizes = np.asarray([len(c) for c in components_cut], dtype=float)

    if len(sizes) == 1:
        return max(float(sizes[0]), 1.0)

    q50 = float(np.percentile(sizes, 50))
    main_sizes = sizes[sizes >= q50]

    if len(main_sizes) == 0:
        main_sizes = sizes

    ref_size = float(np.median(main_sizes))
    return max(ref_size, 1.0)


def compute_raw_component_center_caps(components_raw, ref_size):
    """
    根据 raw MNG 连通分量大小，自动估计每个 raw component 最多允许几个中心。

    例子：
    - 普通 diamond 大小约等于 ref_size -> cap = 1
    - 如果顶部三个 diamond 被 raw MNG 连成一个大 component，
      其大小约为 3 * ref_size -> cap = 3
    - 如果 M_cc 把一个普通 diamond 切成很多碎片，
      raw component 大小仍接近 ref_size -> cap 仍然是 1
    """
    caps = {}

    for cid, comp in enumerate(components_raw):
        size = len(comp)

        cap = int(round(size / (ref_size + EPS)))
        cap = max(1, cap)

        caps[cid] = cap

    return caps






def compute_component_delta(data, rho, component_ids):
    """
    不再计算依赖向量，只计算连通分量内的相对距离 delta。

    对每个点 i：
    1. 只在同一连通分量内寻找更高密度点；
    2. 如果 i 是本连通分量内最高密度点，则 delta 取该分量内最大距离；
    3. 不允许跨连通分量寻找依赖点。

    这样可以避免不同簇的重要高密度点互相连接。
    """
    n = len(data)
    delta = np.zeros(n, dtype=float)

    for i in range(n):
        comp_i = component_ids[i]

        same_comp = np.where(component_ids == comp_i)[0]
        higher = [j for j in same_comp if rho[j] > rho[i]]

        if len(higher) == 0:
            # 本连通分量内最高密度点
            if len(same_comp) <= 1:
                delta[i] = 0.0
            else:
                dists = np.linalg.norm(data[same_comp] - data[i], axis=1)
                delta[i] = float(np.max(dists))
        else:
            dists = [euclidean_distance(data, i, int(j)) for j in higher]
            delta[i] = float(np.min(dists))

    return minmax_normalize(delta)


def gravity_distance(data, sigma, i, j):
    """D^G_ij = d^2 / sigma。若无结构相似性，则退化为较大距离。"""
    dij = euclidean_distance(data, i, j)
    sij = get_sigma(sigma, i, j)
    return (dij * dij + EPS) / (sij + EPS)





def compute_structural_reliability(M, S):
    """R_i = |S_i| / |M_i|。"""
    n = len(M)
    R = np.zeros(n, dtype=float)
    for i in range(n):
        R[i] = len(S[i]) / (len(M[i]) + EPS)
    return R


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
    beta=1.5,
    alpha_min=1.3,
    alpha_max=2.5
):
    """
    点 i 到簇原型 P 的改进引力。

    如果 point_i -> j 在 G_norm 中存在，优先使用图上引力；
    如果不存在，则即时计算改进引力。

    同时加入原型距离权重，避免远距离高密度原型点吸走边界点。
    """
    prototype = list(prototype)

    if len(prototype) == 0:
        return 0.0

    # 点到原型的距离尺度，用中位数更稳健
    dists = []
    for j in prototype:
        j = int(j)
        dists.append(euclidean_distance(data, point_i, j))

    dists = np.asarray(dists, dtype=float)

    if np.any(dists > EPS):
        proto_scale = float(np.median(dists[dists > EPS]))
    else:
        proto_scale = 1.0

    force = 0.0

    for j, dij in zip(prototype, dists):
        j = int(j)

        gij = G_norm.get((point_i, j), None)

        # 如果不是 MNG 图中的直接边，则即时计算改进引力
        if gij is None:
            if reliability is not None and local_scale is not None and sigma is not None:
                sij = get_sigma(sigma, point_i, j)

                if sij <= EPS:
                    gij = 0.0
                else:
                    pair_scale = math.sqrt(local_scale[point_i] * local_scale[j]) + EPS
                    d_scaled = dij / pair_scale

                    alpha_ij = adaptive_alpha(
                        point_i,
                        j,
                        reliability,
                        alpha_min=alpha_min,
                        alpha_max=alpha_max
                    )

                    wd = math.exp(-(d_scaled ** 2))

                    # 局部化 β（同步）
                    rel_ij = min(reliability[point_i], reliability[j])
                    beta_local = beta + 0.3 + (1.0 - rel_ij) * 0.7

                    gij = (sij ** beta_local) * wd / ((d_scaled + EPS) ** alpha_ij)

            else:
                # 兜底：旧版引力
                sij = structure_similarity_from_gamma(gamma, point_i, j)
                gij = sij / (dij * dij + EPS)

        # 新增：点到原型的距离权重
        proto_wd = math.exp(-(dij * dij) / (proto_scale * proto_scale + EPS))

        # 质量项仍然使用密度和逆邻居
        mass = rho[j] * (1.0 + nb[j])

        # 可靠性高的原型点略增强
        if reliability is not None:
            mass = mass * (1.0 + reliability[j])

        force += mass * gij * proto_wd

    return force / (math.sqrt(len(prototype)) + EPS)


def compute_correction_factor(
    data,
    gamma,
    rho,
    nb,
    G_norm,
    prototypes,
    reliability=None,
    local_scale=None,
    sigma=None
):
    """
    校正因子：phi_i = 1 / (1 + max_k F(i,P_k))。
    已经被已有簇强吸引的点，中心分数自动降低。
    """
    n = len(data)
    phi = np.ones(n, dtype=float)
    if len(prototypes) == 0:
        return phi

    for i in range(n):
        max_force = 0.0
        for P in prototypes:
            f = prototype_force_to_point(
                data,
                gamma,
                rho,
                nb,
                G_norm,
                i,
                P,
                reliability=reliability,
                local_scale=local_scale,
                sigma=sigma
            )
            if f > max_force:
                max_force = f
        phi[i] = 1.0 / (1.0 + max_force)
    return phi


def compute_center_score(rho, delta, nb, R, phi):
    """
    Score_i = rho_i * delta_i * (1+nb_i) * (1+R_i) * phi_i
    为避免 nb 量纲过强，内部使用归一化 nb。
    """
    nb_n = minmax_normalize(nb)
    score = rho * delta * (1.0 + nb_n) * (1.0 + R) * phi
    return score


def _scaled_edge_distance(data, local_scale, i, j):
    dij = euclidean_distance(data, i, j)
    return dij / (math.sqrt(local_scale[i] * local_scale[j]) + EPS)


def build_component_mng_by_strong_edges(data, M, S, G, local_scale):
    """
    只用于 weak connected components 的强边图。

    注意：
    - M 仍用于 Layer1/Layer2/Layer3 的候选搜索；
    - M_cc 只用于 component_ids；
    - 这样可以避免一两条弱桥边把多个真实簇合成一个 CC。
    """
    n = len(M)

    # 每个点自己的原始引力阈值。这里用 raw G，不用 G_norm，
    # 因为 G_norm 会被行均值归一化，可能放大孤立弱边。
    row_thr = np.full(n, np.inf, dtype=float)

    for i in range(n):
        vals = [
            G.get((i, int(j)), 0.0)
            for j in M[i]
            if G.get((i, int(j)), 0.0) > EPS
        ]
        if len(vals) > 0:
            row_thr[i] = float(np.percentile(vals, 30))

    # 用 SNN 边的 d_scaled 分布自适应估计“正常簇内距离尺度”
    snn_dscaled = []
    for i in range(n):
        for j in S[i]:
            j = int(j)
            if i == j:
                continue
            snn_dscaled.append(_scaled_edge_distance(data, local_scale, i, j))

    if len(snn_dscaled) > 0:
        tau_mutual = float(np.percentile(snn_dscaled, 95))
        tau_strong = float(np.percentile(snn_dscaled, 85))
    else:
        tau_mutual = 1.8
        tau_strong = 1.5

    M_cc = [set() for _ in range(n)]

    for i in range(n):
        for j in M[i]:
            j = int(j)
            if i == j:
                continue

            d_scaled = _scaled_edge_distance(data, local_scale, i, j)

            raw_ij = G.get((i, j), 0.0)
            raw_ji = G.get((j, i), 0.0)

            mutual = (j in S[i])

            # 互邻边也不能无条件保留：跨簇桥边有时也可能互邻
            keep_mutual = (
                mutual
                and d_scaled <= tau_mutual
                and max(raw_ij, raw_ji) >= min(row_thr[i], row_thr[j])
            )

            # 非互邻边必须双向都强，并且距离尺度更严格
            keep_strong = (
                raw_ij >= row_thr[i]
                and raw_ji >= row_thr[j]
                and d_scaled <= tau_strong
            )

            if keep_mutual or keep_strong:
                M_cc[i].add(j)
                M_cc[j].add(i)

    return M_cc


# ============================================================
# 9. 第一层：高置信簇核扩展
# ============================================================

def layer1_core_expansion(
    center,
    labels,
    current_label,
    M,
    S,
    G_norm,
    rho,
    component_ids,
):
    """
    第一层：放宽后的簇核扩展。

    不再依赖 Dep(j)=i。
    从中心出发，在同一连通分量内扩展满足以下任一条件的点：

    1. 与当前 core 存在严格自然邻居关系；
    2. 与当前 core 的最大引力达到该点局部出边引力的较低阈值；
    3. 密度不低于该连通分量的低分位密度。

    目标：
        让第一层形成一个真正的簇核区域，而不是中心附近几个点。
    """
    n = len(M)
    center = int(center)
    comp_id = component_ids[center]

    comp_points = np.where(component_ids == comp_id)[0]
    rho_floor = np.percentile(rho[comp_points], 10)

    # 每个点自己的局部引力阈值，取 30% 分位，避免过严
    local_g_threshold = np.zeros(n, dtype=float)

    for i in range(n):
        vals = [G_norm.get((i, int(j)), 0.0) for j in M[i]]
        vals = [v for v in vals if v > 0]

        if len(vals) == 0:
            local_g_threshold[i] = 0.0
        else:
            local_g_threshold[i] = float(np.percentile(vals, 30))

    core = set([center])
    frontier = set([center])

    max_rounds = 4

    for _ in range(max_rounds):
        new_nodes = set()

        candidate_set = set()
        for u in frontier:
            candidate_set.update(M[u])

        for j in candidate_set:
            j = int(j)

            if component_ids[j] != comp_id:
                continue

            if labels[j] != -1 and labels[j] != current_label:
                continue

            if j in core:
                continue

            # 与 core 的严格自然邻居关系
            strict_hit = False
            for u in core:
                if j in S[u] or u in S[j]:
                    strict_hit = True
                    break

            # 与 core 的最大引力
            g_to_core = 0.0
            for u in core:
                g_to_core = max(
                    g_to_core,
                    G_norm.get((j, u), 0.0),
                    G_norm.get((u, j), 0.0)
                )

            gravity_hit = g_to_core >= local_g_threshold[j]

            density_hit = rho[j] >= rho_floor

            if strict_hit or (gravity_hit and density_hit):
                new_nodes.add(j)

        if len(new_nodes) == 0:
            break

        core.update(new_nodes)
        frontier = new_nodes

    return core


# ============================================================
# 10. 第二层：簇原型扩展
# ============================================================

def layer2_prototype_expansion(
    center,
    core,
    labels,
    current_label,
    data,
    rho,
    M,
    gamma,
    G_norm,
    component_ids,
):
    """
    第二层：放宽后的簇原型扩展。

    改动：
    1. 候选点从一阶邻域扩大到二阶邻域；
    2. 阈值从 mean 改为 25% 分位；
    3. 允许多轮扩展；
    4. 只限制在同一连通分量内，不跨分量。
    """
    center = int(center)
    comp_id = component_ids[center]

    prototype = set(core)
    max_rounds = 5

    for _ in range(max_rounds):
        candidate_set = set()

        # 一阶邻域
        for u in prototype:
            candidate_set.update(M[u])

        # 二阶邻域
        second_order = set()
        for u in candidate_set:
            second_order.update(M[int(u)])

        candidate_set.update(second_order)

        candidate_set = {
            int(v)
            for v in candidate_set
            if component_ids[int(v)] == comp_id
            and v not in prototype
            and (labels[int(v)] == -1 or labels[int(v)] == current_label)
        }

        if len(candidate_set) == 0:
            break

        gamma_proto = set()
        for u in prototype:
            gamma_proto.update(gamma[u])

        scores = {}

        for i in candidate_set:
            g_values = []

            for u in prototype:
                g_values.append(
                    max(
                        G_norm.get((i, u), 0.0),
                        G_norm.get((u, i), 0.0)
                    )
                )

            g_values = np.asarray(g_values, dtype=float)

            g_max = float(np.max(g_values)) if len(g_values) > 0 else 0.0
            g_mean = float(np.mean(g_values)) if len(g_values) > 0 else 0.0

            inter = len(gamma[i] & gamma_proto)
            union = len(gamma[i] | gamma_proto)
            simi = inter / (union + EPS)

            density_ratio = rho[i] / (rho[center] + EPS)
            density_ratio = min(1.0, density_ratio)

            # 第二层得分：引力为主，结构相似为辅
            scores[i] = (
                0.45 * g_max
                + 0.35 * g_mean
                + 0.20 * simi
            ) * (0.5 + 0.5 * density_ratio)

        if len(scores) == 0:
            break

        score_values = np.asarray(list(scores.values()), dtype=float)

        # 原来用 mean 太严格；这里改成 25% 分位
        threshold = float(np.percentile(score_values, 25))

        new_nodes = {i for i, s in scores.items() if s >= threshold}

        if len(new_nodes) == 0:
            break

        old_size = len(prototype)
        prototype.update(new_nodes)

        if len(prototype) == old_size:
            break

    return prototype


# ============================================================
# 11. 第三层：引力边界分配
# ============================================================

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
    """
    第三层：边界吸收层。

    只处理当前原型周围的一阶/二阶候选点；
    不再试图一次吃掉整个连通分量。
    剩余 clean 点交给后续 DPC 分配。
    """
    prototype = set(prototype)

    if len(prototype) == 0:
        return set()

    ref = next(iter(prototype))
    comp_id = component_ids[int(ref)]

    # 1. 候选点：prototype 的一阶 + 二阶邻域
    candidate_set = set()

    for u in prototype:
        candidate_set.update(M[int(u)])

    second_order = set()
    for u in candidate_set:
        second_order.update(M[int(u)])

    candidate_set.update(second_order)

    candidate_set = {
        int(i)
        for i in candidate_set
        if labels[int(i)] == -1
        and component_ids[int(i)] == comp_id
        and int(i) not in prototype
    }

    if len(candidate_set) == 0:
        return set()

    # 2. 只和同一连通分量内已有原型竞争
    same_comp_prototypes = []
    same_comp_proto_indices = []

    for k, P in enumerate(prototypes):
        if len(P) == 0:
            continue

        p_ref = next(iter(P))
        if component_ids[int(p_ref)] == comp_id:
            same_comp_prototypes.append(P)
            same_comp_proto_indices.append(k)

    all_prototypes = same_comp_prototypes + [set(prototype)]
    current_index = len(all_prototypes) - 1

    # 3. 当前原型对候选点的吸引力
    current_forces = []

    for i in candidate_set:
        f = prototype_force_to_point(
            data,
            gamma,
            rho,
            nb,
            G_norm,
            i,
            prototype,
            reliability=reliability,
            local_scale=local_scale,
            sigma=sigma
        )

        # 如果图引力为 0，用最近原型距离兜底，避免第三层完全推不动
        if f <= EPS:
            dmin = min(euclidean_distance(data, i, int(p)) for p in prototype)
            f = 1.0 / (dmin + EPS)

        current_forces.append(f)

    current_forces = np.asarray(current_forces, dtype=float)

    positive_forces = current_forces[current_forces > EPS]

    if len(positive_forces) == 0:
        return set()

    # 原来 20% 还是偏严，这里改成 5%
    force_threshold = float(np.percentile(positive_forces, 5))

    layer3 = set()

    for i in candidate_set:
        forces = []

        for P in all_prototypes:
            f = prototype_force_to_point(
                data,
                gamma,
                rho,
                nb,
                G_norm,
                i,
                P,
                reliability=reliability,
                local_scale=local_scale,
                sigma=sigma
            )

            if f <= EPS:
                dmin = min(euclidean_distance(data, i, int(p)) for p in P)
                f = 1.0 / (dmin + EPS)

            forces.append(f)

        best_index = int(np.argmax(forces))

        if best_index != current_index:
            continue

        if forces[current_index] >= force_threshold:
            layer3.add(i)

    return layer3


# ============================================================
# 12. clean 数据上的 MNG-DPC 主聚类
# ============================================================
def dpc_assign_remaining_clean_points_v94(
    data,
    labels,
    rho,
    component_ids_cut,
    component_ids_raw,
    raw_parent=None
):
    labels = labels.copy()
    n = len(data)

    if raw_parent is None:
        raw_parent = {}

    order = np.argsort(-rho)

    for i in order:
        i = int(i)

        if labels[i] != -1:
            continue

        cut_i = int(component_ids_cut[i])
        raw_i = int(component_ids_raw[i])
        parent_i = int(raw_parent.get(raw_i, raw_i))

        # 1. 优先同 cut component 内继承
        same_cut = np.where(component_ids_cut == cut_i)[0]

        candidates = [
            int(j)
            for j in same_cut
            if rho[int(j)] > rho[i] and labels[int(j)] != -1
        ]

        if len(candidates) > 0:
            nearest = min(
                candidates,
                key=lambda j: euclidean_distance(data, i, j)
            )
            labels[i] = labels[nearest]
            continue

        same_cut_labeled = [
            int(j)
            for j in same_cut
            if labels[int(j)] != -1
        ]

        if len(same_cut_labeled) > 0:
            nearest = min(
                same_cut_labeled,
                key=lambda j: euclidean_distance(data, i, j)
            )
            labels[i] = labels[nearest]
            continue

        # 2. 再到同 raw component 继承
        same_raw = np.where(component_ids_raw == raw_i)[0]

        same_raw_labeled = [
            int(j)
            for j in same_raw
            if labels[int(j)] != -1
        ]

        if len(same_raw_labeled) > 0:
            nearest = min(
                same_raw_labeled,
                key=lambda j: euclidean_distance(data, i, j)
            )
            labels[i] = labels[nearest]
            continue

        # 3. 如果该 raw component 被吸收，则到 parent raw component 继承
        same_parent = np.where(component_ids_raw == parent_i)[0]

        same_parent_labeled = [
            int(j)
            for j in same_parent
            if labels[int(j)] != -1
        ]

        if len(same_parent_labeled) > 0:
            nearest = min(
                same_parent_labeled,
                key=lambda j: euclidean_distance(data, i, j)
            )
            labels[i] = labels[nearest]
            continue

        # 4. 最后全局最近兜底
        global_labeled = np.where(labels != -1)[0]

        if len(global_labeled) > 0:
            nearest = min(
                global_labeled,
                key=lambda j: euclidean_distance(data, i, int(j))
            )
            labels[i] = labels[nearest]
        else:
            labels[i] = 0

    return labels

def cluster_clean_data(clean_data):
    """
    在 clean 数据上执行 MNG-DPC 主过程。

    返回：
        labels: clean 数据标签，0,1,2,...；未分配点会在函数末尾分配
        prototypes: 每个簇的原型点集合，索引是 clean_data 内部索引
        info: 中间信息
    """
    n = len(clean_data)
    if n == 0:
        return np.asarray([], dtype=int), [], {}

    r, NN, RNN, nb, tree = natural_neighbor_search(clean_data)
    M, S, gamma = build_mng(NN, RNN)

    # 新增：可靠性与局部尺度
    reliability = compute_reliability_from_nb(nb)
    local_scale = compute_local_scale_from_gamma(clean_data, gamma)

    # 改进结构相似性与引力
    sigma = compute_structure_similarity(M, gamma, reliability)
    G, G_norm, mean_G = compute_gravity(
        clean_data,
        M,
        sigma,
        reliability,
        local_scale
    )
    rho, rho_detail = compute_composite_density(clean_data, NN, RNN, M, S, G)
    # 1. raw component：用于中心预算
    component_ids_raw, components_raw = weak_connected_components(M)

    # 2. cut component：用于 delta、扩展约束、DPC 剩余分配
    M_cc = build_component_mng_by_strong_edges(
        clean_data,
        M,
        S,
        G,
        local_scale
    )

    component_ids, components = weak_connected_components(M_cc)

    # 3. 中心预算
    ref_component_size = robust_reference_component_size(components)
    raw_center_caps, raw_cap_info = compute_raw_component_center_caps_v94(
        data=clean_data,
        components_raw=components_raw,
        components_cut=components,
        component_ids_raw=component_ids_raw,
        ref_size=ref_component_size,
        local_scale=local_scale
    )

    max_center_num = int(sum(raw_center_caps.values()))

    # 4. delta 仍然用 M_cc component
    delta = compute_component_delta(clean_data, rho, component_ids)

    R = compute_structural_reliability(M, S)

    base_phi = np.ones(n, dtype=float)
    base_score = compute_center_score(rho, delta, nb, R, base_phi)

    # 5. 先一次性选出中心
    # 5. 先一次性选出中心
    centers, center_records = select_centers_by_raw_cut_budget(
        base_score=base_score,
        components_cut=components,
        component_ids_raw=component_ids_raw,
        components_raw=components_raw,
        raw_center_caps=raw_center_caps
    )

    # 统计每个 raw component 实际拿到几个中心
    raw_center_count = {
        cid: sum(1 for r in center_records if r["raw_cid"] == cid)
        for cid in range(len(components_raw))
    }

    # 只考虑 cap > 0 的 raw component
    active_raw_components = [
        cid for cid, cap in raw_center_caps.items()
        if cap > 0
    ]

    # 有预算但一个中心都没拿到的
    uncovered_raw_components = [
        cid
        for cid in active_raw_components
        if raw_center_count[cid] == 0
    ]

    # 中心数还没用满预算的
    expandable_raw_components = [
        cid
        for cid in active_raw_components
        if raw_center_count[cid] < raw_center_caps[cid]
    ]

    raw_to_cut = build_raw_to_cut_component_map(
        component_ids_raw,
        components
    )

    # 记录已选过的 cut component
    used_cut_ids = set(r["cut_cid"] for r in center_records)

    # 第二轮：先补 uncovered（完全没有中心的可疑簇）
    for raw_cid in uncovered_raw_components:
        cut_ids = sorted(
            raw_to_cut.get(raw_cid, []),
            key=lambda cid: len(components[cid]),
            reverse=True
        )
        for cut_cid in cut_ids:
            idx = np.asarray(list(components[cut_cid]), dtype=int)
            if len(idx) == 0:
                continue
            local_best = int(idx[np.argmax(base_score[idx])])
            center_records.append({
                "center": local_best,
                "raw_cid": raw_cid,
                "cut_cid": cut_cid,
                "cut_size": len(idx),
                "score": float(base_score[local_best]),
            })
            used_cut_ids.add(cut_cid)
            raw_center_count[raw_cid] += 1
            break

    # 第三轮：再补 expandable（预算没用完的）
    for raw_cid in expandable_raw_components:
        remain = raw_center_caps[raw_cid] - raw_center_count[raw_cid]
        if remain <= 0:
            continue
        cut_ids = sorted(
            raw_to_cut.get(raw_cid, []),
            key=lambda cid: len(components[cid]),
            reverse=True
        )
        for cut_cid in cut_ids:
            if cut_cid in used_cut_ids:
                continue
            idx = np.asarray(list(components[cut_cid]), dtype=int)
            if len(idx) == 0:
                continue
            local_best = int(idx[np.argmax(base_score[idx])])
            center_records.append({
                "center": local_best,
                "raw_cid": raw_cid,
                "cut_cid": cut_cid,
                "cut_size": len(idx),
                "score": float(base_score[local_best]),
            })
            used_cut_ids.add(cut_cid)
            raw_center_count[raw_cid] += 1
            remain -= 1
            if remain <= 0:
                break

    centers = [r["center"] for r in center_records]


    labels = np.full(n, -1, dtype=int)
    prototypes = []
    clusters = []
    history = []

    # 6. 先把所有中心预标记，防止前面的中心吞掉后面的中心

    # 6. 先把所有中心预标记，防止前面的中心吞掉后面的中心
    for lab, c in enumerate(centers):
        labels[int(c)] = lab

    # 7. 再围绕这些固定中心做三层扩展
    for current_label, best_center in enumerate(centers):
        best_center = int(best_center)

        core = layer1_core_expansion(
            best_center,
            labels,
            current_label,
            M,
            S,
            G_norm,
            rho,
            component_ids
        )

        if len(core) == 0:
            core = set([best_center])

        prototype = layer2_prototype_expansion(
            best_center,
            core,
            labels,
            current_label,
            clean_data,
            rho,
            M,
            gamma,
            G_norm,
            component_ids
        )

        cluster_points = set(prototype)

        for p in cluster_points:
            labels[int(p)] = current_label

        layer3 = layer3_gravity_boundary_assignment(
            prototype,
            labels,
            current_label,
            prototypes,
            clean_data,
            gamma,
            rho,
            nb,
            M,
            G_norm,
            component_ids,
            reliability=reliability,
            local_scale=local_scale,
            sigma=sigma
        )

        for p in layer3:
            labels[int(p)] = current_label

        cluster_points.update(layer3)

        prototypes.append(set(prototype))
        clusters.append(set(cluster_points))

        history.append({
            "center": int(best_center),
            "core": set(core),
            "prototype": set(prototype),
            "layer3": set(layer3),
            "cluster_points": set(cluster_points),
            "center_record": center_records[current_label],
            "labels_after": labels.copy(),
        })

    labels = dpc_assign_remaining_clean_points_v94(
        data=clean_data,
        labels=labels,
        rho=rho,
        component_ids_cut=component_ids,
        component_ids_raw=component_ids_raw,
        raw_parent=raw_cap_info["raw_parent"]
    )

    info = {
        "r": r,
        "NN": NN,
        "RNN": RNN,
        "nb": nb,
        "M": M,
        "S": S,
        "gamma": gamma,
        "sigma": sigma,
        "G_norm": G_norm,
        "mean_G": mean_G,
        "rho": rho,
        "rho_detail": rho_detail,
        "component_ids": component_ids,
        "components": components,
        "delta": delta,
        "R": R,
        "base_score": base_score,

        "history": history,

        # 新增
        "reliability": reliability,
        "local_scale": local_scale,
        "M_cc": M_cc,
        "component_ids_raw": component_ids_raw,
        "components_raw": components_raw,
        "ref_component_size": ref_component_size,
        "raw_center_caps": raw_center_caps,
        "ref_component_size": ref_component_size,
        "raw_center_caps": raw_center_caps,
        "max_center_num": len(centers),
        "raw_center_count": {
            cid: sum(1 for r in center_records if r["raw_cid"] == cid)
            for cid in range(len(components_raw))
        },

        "raw_cap_info": raw_cap_info,


    }
    return labels.astype(int), prototypes, info


# ============================================================
# 13. 退化方案：MNG 连通分量
# ============================================================

def fallback_connected_components(M):
    n = len(M)
    labels = np.full(n, -1, dtype=int)
    prototypes = []
    clusters = []
    history = []
    current_label = 0

    for start in range(n):
        if labels[start] != -1:
            continue
        comp = set()
        q = deque([start])
        labels[start] = current_label
        while q:
            i = q.popleft()
            comp.add(i)
            for j in M[i]:
                if labels[j] == -1:
                    labels[j] = current_label
                    q.append(j)
        prototypes.append(set(comp))
        clusters.append(set(comp))
        current_label += 1

    return labels, prototypes, clusters


# ============================================================
# 14. 原始数据低密度点重分配
# ============================================================

def compute_original_graph_info(data):
    """原始数据上计算 MNG、密度、引力等，用于低密度点重分配。"""
    r, NN, RNN, nb, tree = natural_neighbor_search(data)
    M, S, gamma = build_mng(NN, RNN)

    # 新增：可靠性与局部尺度
    reliability = compute_reliability_from_nb(nb)
    local_scale = compute_local_scale_from_gamma(data, gamma)

    # 改进结构相似性与引力
    sigma = compute_structure_similarity(M, gamma, reliability)
    G, G_norm, mean_G = compute_gravity(
        data,
        M,
        sigma,
        reliability,
        local_scale
    )
    rho, rho_detail = compute_composite_density(data, NN, RNN, M, S, G)
    Q = compute_low_confidence_score(rho, nb, mean_G)
    return {
        "r": r,
        "NN": NN,
        "RNN": RNN,
        "nb": nb,
        "M": M,
        "S": S,
        "gamma": gamma,
        "sigma": sigma,
        "G_norm": G_norm,
        "mean_G": mean_G,
        "rho": rho,
        "rho_detail": rho_detail,
        "Q": Q,
        "tree": tree,

        # 新增
        "reliability": reliability,
        "local_scale": local_scale,
    }


def redistribute_low_density_points(data, labels, low_indices, prototypes_original, original_info):
    """
    对低密度点做点对簇原型引力重分配。
    """
    gamma = original_info["gamma"]
    rho = original_info["rho"]
    nb = original_info["nb"]
    G_norm = original_info["G_norm"]

    reliability = original_info.get("reliability", None)
    local_scale = original_info.get("local_scale", None)
    sigma = original_info.get("sigma", None)

    if len(prototypes_original) == 0:
        return labels

    assigned_indices = np.where(labels != -1)[0]

    for i in low_indices:
        i = int(i)
        forces = [
            prototype_force_to_point(
                data,
                gamma,
                rho,
                nb,
                G_norm,
                i,
                P,
                reliability=reliability,
                local_scale=local_scale,
                sigma=sigma
            )
            for P in prototypes_original
        ]
        if len(forces) > 0 and max(forces) > 0:
            labels[i] = int(np.argmax(forces))
        else:
            # 若结构引力为 0，则退化到最近已分配点。
            if len(assigned_indices) == 0:
                labels[i] = 0
            else:
                nearest = min(assigned_indices, key=lambda j: euclidean_distance(data, i, int(j)))
                labels[i] = labels[nearest]

    return labels


# ============================================================
# 15. MNG-DPC 总入口
# ============================================================

def mng_dpc(data):
    """
    MNG-DPC 总流程。

    返回：
        labels: 原始数据标签
        info: 中间信息
    """
    data = np.asarray(data, dtype=float)
    n = len(data)
    if n == 0:
        return np.asarray([], dtype=int), {}

    # 1. 原始图信息
    original_info = compute_original_graph_info(data)

    # 2. 自动识别低密度点：DSNGCAP + nb自适应 双方案并集
    low_ds, clean_ds, alpha_ds, denoise_ds = detect_low_density_points_by_dsngcap_density(
        data=data,
        tree=original_info["tree"],
        r=original_info["r"],
        nb=original_info["nb"]
    )
    low_nb, clean_nb, alpha_nb, denoise_nb = detect_low_density_points_by_nb_adaptive(
        nb=original_info["nb"],
        r=original_info["r"]
    )

    # 取并集：两个方案任一认定为低密度的点都去掉
    low_indices = np.union1d(low_ds, low_nb).astype(int)
    clean_indices = np.setdiff1d(np.arange(n, dtype=int), low_indices).astype(int)
    alpha = int(round(len(low_indices) / n * 100))

    denoise_info = {
        "dsngcap": denoise_ds,
        "nb_adaptive": denoise_nb,
        "union_low_count": len(low_indices),
        "dsngcap_only": len(np.setdiff1d(low_indices, low_nb)),
        "nb_only": len(np.setdiff1d(low_nb, low_ds)),
        "both": len(np.intersect1d(low_ds, low_nb)),
    }

    clean_data = data[clean_indices]

    # 3. clean 数据上重新建图并执行三层 MNG-DPC
    clean_labels, clean_prototypes, clean_info = cluster_clean_data(clean_data)

    # 4. clean 标签映射回原始数据
    labels = np.full(n, -1, dtype=int)
    for local_i, original_i in enumerate(clean_indices):
        labels[int(original_i)] = int(clean_labels[local_i])

    # 5. 原型点索引映射回原始数据
    prototypes_original = []
    for P in clean_prototypes:
        prototypes_original.append(set(int(clean_indices[p]) for p in P))

    # 6. 低密度点引力重分配
    labels = redistribute_low_density_points(data, labels, low_indices, prototypes_original, original_info)

    info = {
        "alpha": alpha,
        "low_indices": low_indices,
        "clean_indices": clean_indices,
        "original_info": original_info,
        "clean_info": clean_info,
        "clean_labels": clean_labels,
        "clean_prototypes": clean_prototypes,
        "prototypes_original": prototypes_original,
        "num_clusters": len(set(labels.tolist())) if len(labels) > 0 else 0,
        "final_labels": labels.copy(),
        "denoise_info": denoise_info,
    }
    return labels.astype(int), info


# ============================================================
# 16. 可视化与简单输出
# ============================================================

def save_and_close(save_path, dpi=160, show=False):
    """
    快速保存图片，并立即关闭，避免 PyCharm / matplotlib 长时间占用内存。
    """
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

def plot_fast_final_result(data, labels, title, save_path=None, show=False):
    """
    最省时的最终聚类结果图：
    1. 不画边；
    2. 不画中心；
    3. 不画额外文字；
    4. 只画不同类不同颜色。
    """
    data_2d = project_to_2d(data)
    labels = np.asarray(labels, dtype=int)

    plt.figure(figsize=(7, 6))

    unique_labels = sorted(set(labels.tolist()))

    for lab in unique_labels:
        idx = np.where(labels == lab)[0]
        plt.scatter(
            data_2d[idx, 0],
            data_2d[idx, 1],
            s=6,
            alpha=0.85,
            label=f"C{lab}"
        )

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")

    if len(unique_labels) <= 20:
        plt.legend(markerscale=2, fontsize=8)

    save_and_close(save_path, dpi=140, show=show)

def plot_final_result_with_centers(data, labels, info, title, save_path=None, show=False):
    """
    最终结果 + 所有中心编号。
    不再每个中心单独画一张图。
    """
    data_2d = project_to_2d(data)
    labels = np.asarray(labels, dtype=int)

    clean_indices = info["clean_indices"]
    history = info["clean_info"].get("history", [])

    # clean_data 内部中心索引 -> 原始数据索引
    center_original_indices = []
    for t, step in enumerate(history, start=1):
        c_local = int(step["center"])
        c_original = int(clean_indices[c_local])
        center_original_indices.append((t, c_original))

    plt.figure(figsize=(7, 6))

    unique_labels = sorted(set(labels.tolist()))

    for lab in unique_labels:
        idx = np.where(labels == lab)[0]
        plt.scatter(
            data_2d[idx, 0],
            data_2d[idx, 1],
            s=6,
            alpha=0.80,
            label=f"C{lab}"
        )

    # 统一画中心，并标序号
    for t, c in center_original_indices:
        x, y = data_2d[c, 0], data_2d[c, 1]

        plt.scatter(
            x,
            y,
            s=90,
            c="red",
            marker="*",
            edgecolors="black",
            linewidths=0.7,
            zorder=5
        )

        plt.text(
            x,
            y,
            str(t),
            fontsize=9,
            color="black",
            ha="left",
            va="bottom",
            zorder=6
        )

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")

    if len(unique_labels) <= 20:
        plt.legend(markerscale=2, fontsize=8)

    save_and_close(save_path, dpi=160, show=show)

def plot_fast_directed_mng_graph(
    data,
    M,
    NN,
    RNN,
    title,
    save_path=None,
    show=False,
    max_edges=120000
):
    """
    快速 MNG 单向图，邻居(NN)和逆邻居(RNN)分色显示。

    蓝色箭头 = NN（点 i 主动选择的自然邻居）
    红色箭头 = RNN（把 i 当邻居的点，即反向边）
    绿色箭头 = SNN（互为自然邻居，交集部分）

    箭头比原版更粗、更大、更明显。
    """
    data_2d = project_to_2d(data)

    X = data_2d[:, 0]
    Y = data_2d[:, 1]

    # 分别收集 NN、RNN-only、SNN（交集）的边
    nn_x0, nn_y0, nn_dx, nn_dy = [], [], [], []
    rnn_x0, rnn_y0, rnn_dx, rnn_dy = [], [], [], []
    snn_x0, snn_y0, snn_dx, snn_dy = [], [], [], []

    for i in range(len(M)):
        nn_i = set(int(v) for v in NN[i])
        rnn_i = set(int(v) for v in RNN[i])
        snn_i = nn_i & rnn_i

        # NN-only（i 主动选的邻居，对方没选 i）
        for j in nn_i - snn_i:
            nn_x0.append(X[i]); nn_y0.append(Y[i])
            nn_dx.append(X[j] - X[i]); nn_dy.append(Y[j] - Y[i])

        # RNN-only（对方选了 i，但 i 没选对方）
        for j in rnn_i - snn_i:
            rnn_x0.append(X[i]); rnn_y0.append(Y[i])
            rnn_dx.append(X[j] - X[i]); rnn_dy.append(Y[j] - Y[i])

        # SNN（互为自然邻居）
        for j in snn_i:
            snn_x0.append(X[i]); snn_y0.append(Y[i])
            snn_dx.append(X[j] - X[i]); snn_dy.append(Y[j] - Y[i])

    total_edges = len(nn_x0) + len(rnn_x0) + len(snn_x0)

    if total_edges == 0:
        plt.figure(figsize=(10, 9))
        plt.scatter(X, Y, s=5, alpha=0.85)
        plt.title(title + " | no edges")
        plt.xlabel("x"); plt.ylabel("y")
        save_and_close(save_path, dpi=300, show=show)
        return

    # 大图限流：按比例从三类边中抽样
    if max_edges is not None and total_edges > max_edges:
        rng = np.random.default_rng(42)

        def sample_edges(x0, y0, dx, dy, n_total, n_sample):
            arr = list(zip(x0, y0, dx, dy))
            if len(arr) == 0:
                return [], [], [], []
            k = max(0, int(len(arr) / n_total * n_sample))
            k = min(k, len(arr))
            idx = rng.choice(len(arr), size=k, replace=False)
            sampled = [arr[i] for i in idx]
            return [s[0] for s in sampled], [s[1] for s in sampled], [s[2] for s in sampled], [s[3] for s in sampled]

        nn_x0, nn_y0, nn_dx, nn_dy = sample_edges(nn_x0, nn_y0, nn_dx, nn_dy, total_edges, max_edges)
        rnn_x0, rnn_y0, rnn_dx, rnn_dy = sample_edges(rnn_x0, rnn_y0, rnn_dx, rnn_dy, total_edges, max_edges)
        snn_x0, snn_y0, snn_dx, snn_dy = sample_edges(snn_x0, snn_y0, snn_dx, snn_dy, total_edges, max_edges)

    plt.figure(figsize=(10, 9))
    ax = plt.gca()


    # 先画点
    ax.scatter(X, Y, s=4, c="tab:blue", alpha=0.70, zorder=2)

    # 公用箭头参数
    quiver_kw = dict(
        angles="xy", scale_units="xy", scale=1,
        width=0.0015, headwidth=3.5, headlength=4.5, headaxislength=3.8,
        zorder=1
    )

    # NN-only：蓝色
    if len(nn_x0) > 0:
        ax.quiver(nn_x0, nn_y0, nn_dx, nn_dy,
                  color="#1f77b4", alpha=0.40, label=f"NN ({len(nn_x0)} edges)", **quiver_kw)

    # RNN-only：红色
    if len(rnn_x0) > 0:
        ax.quiver(rnn_x0, rnn_y0, rnn_dx, rnn_dy,
                  color="#d62728", alpha=0.40, label=f"RNN ({len(rnn_x0)} edges)", **quiver_kw)

    # SNN（互为邻居）：绿色
    if len(snn_x0) > 0:
        ax.quiver(snn_x0, snn_y0, snn_dx, snn_dy,
                  color="#2ca02c", alpha=0.50, label=f"SNN ({len(snn_x0)} edges)", **quiver_kw)

    if max_edges is not None and total_edges > max_edges:
        title = f"{title} | sampled ~{max_edges}/{total_edges} directed edges"
    else:
        title = f"{title} | {total_edges} directed edges"

    ax.set_title(title)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.legend(markerscale=3, fontsize=8, loc="upper right")
    save_and_close(save_path, dpi=300, show=show)





def project_to_2d(data):
    data = np.asarray(data, dtype=float)
    if data.shape[1] == 1:
        return np.hstack([data, np.zeros((len(data), 1))])
    if data.shape[1] == 2:
        return data[:, :2]
    pca = PCA(n_components=2)
    return pca.fit_transform(data)


def visualize_lightweight(data, labels, info, data_path=None, out_dir=None, show=False):
    """
    轻量版可视化。

    只生成三类图：
    1. 最终聚类结果图，不同类不同颜色；
    2. 最终聚类结果 + 中心编号；
    3. clean MNG 单向图。

    不画：
    - 每个中心单独图；
    - 每轮 Step 11 / Step 12；
    - weak connected components；
    - 密度图；
    - 低置信度图。
    """
    if out_dir is None:
        stem = os.path.splitext(os.path.basename(data_path if data_path is not None else "dataset"))[0]
        out_dir = os.path.join(os.path.dirname(__file__), "实验数据存放地址", f"{stem}_mng_dpc_visuals_light")

    os.makedirs(out_dir, exist_ok=True)

    clean_indices = info["clean_indices"]
    clean_data = data[clean_indices]
    clean_info = info["clean_info"]

    # 1. 最开始先画最终聚类结果图
    plot_fast_final_result(
        data=data,
        labels=labels,
        title="Final clustering result",
        save_path=os.path.join(out_dir, "00_final_result_fast.png"),
        show=show
    )

    # 2. 最终结果 + 所有中心编号
    plot_final_result_with_centers(
        data=data,
        labels=labels,
        info=info,
        title="Final clustering result with center indices",
        save_path=os.path.join(out_dir, "01_final_result_with_centers.png"),
        show=show
    )

    # 3. clean MNG 单向图（蓝=NN, 红=RNN, 绿=互为邻居）
    plot_fast_directed_mng_graph(
        data=clean_data,
        M=clean_info["M"],
        NN=clean_info["NN"],
        RNN=clean_info["RNN"],
        title="Clean MNG directed graph (blue=NN, red=RNN, green=SNN)",
        save_path=os.path.join(out_dir, "02_clean_mng_directed_graph.png"),
        show=show,
        max_edges=120000
    )

    print(f"轻量版可视化图片已保存到: {out_dir}")
    return out_dir


def print_summary(labels, info, elapsed):
    counter = Counter(labels.tolist())
    print('=' * 60)
    print('MNG-DPC finished')
    print(f"低密度/低置信点比例 alpha: {info['alpha']}%")
    print(f"低密度/低置信点数量: {len(info['low_indices'])}")
    print(f"最终簇个数: {len(counter)}")
    print(f"各簇样本数: {dict(counter)}")
    print(f"耗时: {elapsed:.2f} 秒")
    print('=' * 60)


def print_v94_center_budget_diagnostic(clean_info):
    components = clean_info["components"]
    components_raw = clean_info["components_raw"]
    raw_center_caps = clean_info["raw_center_caps"]
    raw_center_count = clean_info["raw_center_count"]

    print("=" * 70)
    print("v9.4 center budget diagnostic")
    print(f"M_cc component count: {len(components)}")
    print(f"Raw M component count: {len(components_raw)}")
    print(f"ref_component_size: {clean_info['ref_component_size']:.2f}")
    print(f"max_center_num: {clean_info['max_center_num']}")

    print("M_cc component sizes:")
    print(sorted([len(c) for c in components]))

    print("Raw component sizes and caps:")
    for cid, comp in enumerate(components_raw):
        print(
            f"raw CC{cid}: size={len(comp)}, "
            f"cap={raw_center_caps[cid]}, "
            f"actual_centers={raw_center_count[cid]}"
        )

    print("=" * 70)
# ============================================================
# 17. 数据集路径放最后
# ============================================================

if __name__ == "__main__":
    start_time = time.time()

    data_path = "Synthetic/diamond9.csv"

    data = load_data(data_path)
    labels, info = mng_dpc(data)

    print_v94_center_budget_diagnostic(info["clean_info"])

    elapsed = time.time() - start_time
    print_summary(labels, info, elapsed)

    # 轻量版可视化：
    # 1. 先画最终聚类结果；
    # 2. 再画中心编号图；
    # 3. 再画 clean MNG 单向图；
    # 4. 不画 weak connected components；
    # 5. 不再单独画每个中心的图。
    visualize_lightweight(data, labels, info, data_path=data_path, show=False)
