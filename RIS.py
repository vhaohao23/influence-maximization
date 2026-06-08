from __future__ import annotations

import argparse
import random
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd


DEFAULT_GRAPH_PATH = Path("network/Wiki-Vote.txt")


Adjacency = Dict[int, List[int]]


def load_graph(path: str | Path = DEFAULT_GRAPH_PATH) -> pd.DataFrame:
    """
    Doc danh sach canh co huong tu file text.

    File duoc gia dinh co 2 cot: source target. Cac dong comment bat dau bang '#'
    se duoc bo qua neu co.
    """
    graph = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["source", "target"],
        usecols=[0, 1],
        dtype=np.int64,
        engine="python",
    )
    return graph.drop_duplicates().reset_index(drop=True)


def build_adjacency(G: pd.DataFrame) -> Tuple[List[int], Adjacency, Adjacency]:
    """
    Tao danh sach node, in-neighbor va out-neighbor tu dataframe canh.
    """
    required_columns = {"source", "target"}
    missing_columns = required_columns - set(G.columns)
    if missing_columns:
        raise ValueError(f"G must contain columns {sorted(required_columns)}")

    in_neighbors: Adjacency = defaultdict(list)
    out_neighbors: Adjacency = defaultdict(list)
    nodes: Set[int] = set()

    for source, target in G[["source", "target"]].itertuples(index=False, name=None):
        source = int(source)
        target = int(target)
        nodes.add(source)
        nodes.add(target)
        out_neighbors[source].append(target)
        in_neighbors[target].append(source)

    if not nodes:
        raise ValueError("G must contain at least one edge")

    return list(nodes), dict(in_neighbors), dict(out_neighbors)


def get_rr_set(
    nodes: Sequence[int],
    in_neighbors: Adjacency,
    p: float,
    rng: random.Random,
) -> Set[int]:
    """
    Sinh mot Reverse Reachable Set theo mo hinh Independent Cascade.

    Chon ngau nhien node dich v, sau do di nguoc cac canh u -> v. Moi canh duoc
    chap nhan voi xac suat p. Cac node di nguoc den duoc v tao thanh RR set.
    """
    if not 0 <= p <= 1:
        raise ValueError("p must be in [0, 1]")

    start = rng.choice(nodes)
    rr_set = {start}
    frontier = [start]

    while frontier:
        current = frontier.pop()
        for predecessor in in_neighbors.get(current, []):
            if predecessor not in rr_set and rng.random() < p:
                rr_set.add(predecessor)
                frontier.append(predecessor)

    return rr_set


def _select_seeds_from_rr_sets(
    rr_sets: List[Set[int]],
    k: int,
    nodes: Sequence[int],
) -> List[int]:
    """
    Greedy maximum coverage tren tap RR sets.
    """
    if k <= 0:
        return []

    rr_by_node: Dict[int, Set[int]] = defaultdict(set)
    coverage_count: Counter[int] = Counter()

    for rr_index, rr_set in enumerate(rr_sets):
        for node in rr_set:
            rr_by_node[node].add(rr_index)
            coverage_count[node] += 1

    selected: List[int] = []
    selected_set: Set[int] = set()
    uncovered_rr_sets = set(range(len(rr_sets)))

    for _ in range(min(k, len(nodes))):
        if coverage_count:
            seed, _ = coverage_count.most_common(1)[0]
        else:
            seed = next(node for node in nodes if node not in selected_set)

        selected.append(seed)
        selected_set.add(seed)

        covered_now = rr_by_node.get(seed, set()) & uncovered_rr_sets
        uncovered_rr_sets.difference_update(covered_now)

        for rr_index in covered_now:
            for node in rr_sets[rr_index]:
                if node not in selected_set:
                    coverage_count[node] -= 1
                    if coverage_count[node] <= 0:
                        del coverage_count[node]

        coverage_count.pop(seed, None)

    return selected


def ris(
    G: pd.DataFrame,
    k: int,
    p: float = 0.1,
    mc: int = 1000,
    seed: int | None = None,
) -> Tuple[List[int], List[float]]:
    """
    Reverse Influence Sampling cho bai toan Influence Maximization.

    Args:
        G: DataFrame canh co huong voi columns ['source', 'target'].
        k: So seed can chon.
        p: Xac suat kich hoat cua moi canh theo Independent Cascade.
        mc: So RR sets can sinh. Gia tri lon hon cho ket qua on dinh hon.
        seed: Random seed de tai lap ket qua.

    Returns:
        (seed_set, timelapse), trong do timelapse[i] la thoi gian sau khi chon
        seed thu i + 1.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    if mc <= 0:
        raise ValueError("mc must be positive")

    start_time = time.time()
    rng = random.Random(seed)
    nodes, in_neighbors, _ = build_adjacency(G)

    rr_sets = [get_rr_set(nodes, in_neighbors, p, rng) for _ in range(mc)]
    seeds = _select_seeds_from_rr_sets(rr_sets, k, nodes)
    timelapse = [time.time() - start_time for _ in seeds]

    return sorted(seeds), timelapse


def independent_cascade(
    G: pd.DataFrame,
    S: Sequence[int],
    p: float = 0.1,
    mc: int = 1000,
    seed: int | None = None,
) -> float:
    """
    Uoc luong spread trung binh cua seed set bang Monte Carlo Independent Cascade.
    """
    if not 0 <= p <= 1:
        raise ValueError("p must be in [0, 1]")
    if mc <= 0:
        raise ValueError("mc must be positive")

    rng = random.Random(seed)
    _, _, out_neighbors = build_adjacency(G)
    seed_set = set(map(int, S))
    spreads: List[int] = []

    for _ in range(mc):
        active = set(seed_set)
        frontier = list(seed_set)

        while frontier:
            current = frontier.pop()
            for neighbor in out_neighbors.get(current, []):
                if neighbor not in active and rng.random() < p:
                    active.add(neighbor)
                    frontier.append(neighbor)

        spreads.append(len(active))

    return float(np.mean(spreads))


def IC(
    G: pd.DataFrame,
    S: Sequence[int],
    p: float = 0.1,
    mc: int = 1000,
    seed: int | None = None,
) -> float:
    """
    Alias giu tuong thich voi code cu.
    """
    return independent_cascade(G, S, p=p, mc=mc, seed=seed)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RIS for influence maximization.")
    parser.add_argument("--graph", default=str(DEFAULT_GRAPH_PATH), help="Path to edge list file")
    parser.add_argument("-k", type=int, default=100, help="Number of seeds")
    parser.add_argument("-p", type=float, default=0.1, help="Propagation probability")
    parser.add_argument("--rr", type=int, default=1000, help="Number of RR sets")
    parser.add_argument("--mc", type=int, default=1000, help="IC Monte Carlo runs")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    G = load_graph(args.graph)
    seed_set, timelapse = ris(G, k=args.k, p=args.p, mc=args.rr, seed=args.seed)
    spread = independent_cascade(G, seed_set, p=args.p, mc=args.mc, seed=args.seed + 1)

    print(f"Seeds ({len(seed_set)}): {seed_set}")
    print(f"Estimated spread: {spread:.4f}")
    if timelapse:
        print(f"RIS time: {timelapse[-1]:.4f}s")


if __name__ == "__main__":
    main()
