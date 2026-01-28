limport os
import pandas as pd
from itertools import combinations
from collections import deque

def column_name_to_frozen_coordinate(pdb, csv_path, cutoff=10.0):
    """
    Reads pairwise distances from CSV columns named like 'A_B', builds a graph
    with edges for distances < cutoff, runs BFS to find connectivity, and reports:
      - 'good' + max edge distance (if fully connected), or
      - 'disconnected' + max bridge distance between main component and others.
    """
    df = pd.read_csv(csv_path)
    if len(df) == 0:
        raise ValueError("CSV is empty.")
    if row_index >= len(df):
        raise IndexError(f"row_index {row_index} out of range for CSV with {len(df)} rows.")
    #row = df.iloc[row_index]

    # --- collect chain IDs from columns like 'A_B' ---
    pair_cols = [c for c in df.columns if c.count("_") == 1]
    chains = sorted(set(sum([c.split("_") for c in pair_cols], [])))

    # Optional: if pdb is a real file, use chain IDs from that PDB instead
    if os.path.isfile(pdb):
        try:
            from Bio.PDB import PDBParser
            structure = PDBParser(QUIET=True).get_structure("prot", pdb)
            chains = [ch.id for ch in structure.get_chains()]
        except Exception:
            pass  # fall back to CSV-derived chains if PDB parsing fails

    if not chains:
        raise ValueError("No chain IDs found (from CSV columns or PDB).")

    # --- distances dict keyed by frozenset({'A','B'}) ---
    distances = {}
    for u, v in combinations(chains, 2):
        key, key_r = f"{u}_{v}", f"{v}_{u}"
        val = row[key] if key in row.index else (row[key_r] if key_r in row.index else None)
        if pd.notna(val):
            try:
                distances[frozenset((u, v))] = float(val)
            except (TypeError, ValueError):
                pass

    # --- adjacency: edge if distance < cutoff ---
    adj = {u: set() for u in chains}
    for u, v in combinations(chains, 2):
        d = distances.get(frozenset((u, v)))
        if d is not None and d < cutoff:
            adj[u].add(v); adj[v].add(u)

    # --- BFS main component from first chain ---
    start = chains[0]
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in seen:
                seen.add(w); q.append(w)

    disconnected = [c for c in chains if c not in seen]

    if not disconnected:
        # fully connected under cutoff: report max edge distance among present edges
        edge_dists = []
        for u in chains:
            for v in adj[u]:
                if u < v:
                    d = distances.get(frozenset((u, v)))
                    if d is not None:
                        edge_dists.append(d)
        max_edge = max(edge_dists) if edge_dists else None
        print("good", max_edge)
        return {"status": "good", "max_edge_distance": max_edge, "adjacency": adj, "distances": distances}
    else:
        # disconnected: report max distance among *existing* pairs bridging main->disconnected
        bridge_dists = []
        for u in seen:
            for v in disconnected:
                d = distances.get(frozenset((u, v)))
                if d is not None:
                    bridge_dists.append(d)
        max_bridge = max(bridge_dists) if bridge_dists else None
        print("disconnected", max_bridge)
        return {"status": "disconnected", "max_bridge_distance": max_bridge,
                "main_component": sorted(seen), "disconnected": disconnected,
                "adjacency": adj, "distances": distances}

# Example:

csv = "distance_1AQF_consolidated_score_file.csv"
column_name_to_frozen_coordinate("1AQF", csv)