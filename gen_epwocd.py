"""GEN-EPWOCD community detection (Python wrapper).

Wraps the optimized C++ EP-WOCD implementation that lives next to this file
(``epwocd.cpp`` -> compiled ``epwocd`` binary) behind a NetworkX-friendly API.

The C++ program reads a graph from stdin::

    N NE
    u v        (NE undirected edges, nodes numbered 1..N)

and prints N lines to stdout, line i being the community label (1..k) of node i.
The achieved modularity is printed to stderr.

Public API:
    detect_communities(graph, seed=None) -> (communities, partition, method)
        Drop-in replacement for the notebook's original detect_communities().
    run_labels(graph, seed=None) -> dict[node, int]
        Lower level: returns the raw node -> community-label mapping.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import networkx as nx

_DIR = Path(__file__).resolve().parent
_SRC = _DIR / "epwocd.cpp"
_BIN = _DIR / "epwocd"


def _ensure_built() -> Path:
    """Compile the C++ source if the binary is missing or out of date."""
    if _BIN.exists() and _BIN.stat().st_mtime >= _SRC.stat().st_mtime:
        return _BIN
    if not _SRC.exists():
        raise FileNotFoundError(f"Missing EP-WOCD source: {_SRC}")
    subprocess.run(
        ["g++", "-O2", "-std=c++17", "-o", str(_BIN), str(_SRC)],
        check=True,
    )
    return _BIN


def run_labels(graph: nx.Graph, seed=None) -> dict:
    """Run EP-WOCD on ``graph``; return ``{node: community_label}`` (labels 1..k).

    Node identities are preserved (the graph may use arbitrary hashable nodes);
    they are mapped to consecutive integers only for the C++ process.
    """
    nodes = list(graph.nodes())
    n = len(nodes)
    if n == 0:
        return {}

    node_to_id = {node: i + 1 for i, node in enumerate(nodes)}  # 1-indexed

    lines = [f"{n} {graph.number_of_edges()}"]
    for u, v in graph.edges():
        lines.append(f"{node_to_id[u]} {node_to_id[v]}")
    stdin_data = "\n".join(lines) + "\n"

    binary = _ensure_built()
    args = [str(binary)]
    if seed is not None:
        args.append(str(int(seed)))

    proc = subprocess.run(args, input=stdin_data, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"epwocd failed (exit {proc.returncode}): {proc.stderr.strip()}"
        )

    out = proc.stdout.split()
    if len(out) != n:
        raise RuntimeError(
            f"epwocd returned {len(out)} labels, expected {n}. stderr: {proc.stderr.strip()}"
        )

    return {nodes[i]: int(out[i]) for i in range(n)}


def detect_communities(graph: nx.Graph, seed=None):
    """Drop-in replacement for the notebook's ``detect_communities``.

    Returns ``(communities, partition, method)`` where:
      * ``communities`` -- list[list[node]] sorted by size (desc), nodes sorted
        within each community (matching the original function's output shape);
      * ``partition``   -- dict mapping each node to its community index;
      * ``method``      -- the string ``'gen_epwocd'``.
    """
    labels = run_labels(graph, seed)

    groups: dict = {}
    for node, cid in labels.items():
        groups.setdefault(cid, []).append(node)

    communities = [sorted(comm) for comm in groups.values()]
    communities.sort(key=len, reverse=True)
    partition = {node: cid for cid, comm in enumerate(communities) for node in comm}
    return communities, partition, "gen_epwocd"
