"""Phase I.D: adaptive-network parameter sweep.

Dynamic traditional-network parameter sweep, originally run locally.
Results are saved relative to the current working directory.
The parameter sweep appends results; use a new output file for independent reruns.
"""

import csv
import itertools
import os
import random
from collections import defaultdict

import networkx as nx
import numpy as np
from networkx import grid_graph
from community import community_louvain
from tqdm import tqdm

def mle_tau_fixed_xmin(sizes, s_min):
    tail = np.asarray([s for s in sizes if s >= s_min], dtype=float)
    n = tail.size
    if n == 0:
        return np.nan, 0
    tau = 1.0 + n / np.sum(np.log(tail / s_min))
    return float(tau), int(n)

def ks_distance_fixed_xmin(sizes, s_min, tau):
    tail = np.asarray([s for s in sizes if s >= s_min], dtype=float)
    n = tail.size
    if n == 0 or not np.isfinite(tau) or tau <= 1.0:
        return np.nan
    x = np.sort(tail)
    F_emp = np.arange(1, n + 1) / n
    F_model = 1.0 - (x / s_min) ** (-(tau - 1.0))
    return float(np.max(np.abs(F_emp - F_model)))


def simulate_dynamic_network(N, p, T, A, f, TH):  # f = # of rewirings, TH = threshold of rewiring
    stable_time = max(1, T // f)  # FIX: prevent stable_time=0 when f>T; keeps integer windows

    # --- Build lattice ---
    G = grid_graph(dim=[N, N])
    pos = {(x, y): (x, y) for x, y in G.nodes()}
    nodes = list(G.nodes())  # small speed win
    edges = list(G.edges())

    # --- SW-style rewiring (does NOT count toward persistence) ---
    edges = list(G.edges())
    for edge in edges:
        if random.random() < p:
            u, v = edge
            G.remove_edge(u, v)
            while True:
                new_v = (random.randint(0, N - 1), random.randint(0, N - 1))
                if new_v != u and not G.has_edge(u, new_v):
                    G.add_edge(u, new_v)
                    break

    # --- Initialise info ---
    nx.set_node_attributes(G, {node: random.uniform(0, 1) for node in nodes}, 'info')
    nx.set_edge_attributes(G, 0, 'info_transmitted')  # cumulative flow on edges (within stable window)
    nx.set_edge_attributes(G, 0, 'rewiring_binary')   # 0 if never rewired during Phase I.D, 1 if rewired

    # --- History initialisation ---
    history = []
    branching_history = defaultdict(list)  # per global timestep t
    dynamic_history = []

    edges = list(G.edges())  # refresh after SW rewiring
    total_info = sum(G.nodes[n]['info'] for n in nodes)
    total_info_transmitted = sum(G.edges[e]['info_transmitted'] for e in edges)

    history.append({
        't': 0,
        'avalanche_size': 0,
        'avalanche_duration': None,
        'total_info': total_info,
        'total_info_transmitted': total_info_transmitted,
        'max_wave_size': 0,
        'total_topplings': 0,
        'system_spanning_flag': 0,
        'branching_mean': None,
        'branching_std': None,
        'branching_max': None,
    })

    w = 0  # how many Phase I.D rewires performed
    C = nx.average_clustering(G)
    TR = nx.transitivity(G)
    AS = nx.degree_assortativity_coefficient(G)
    partition = community_louvain.best_partition(G)
    M = community_louvain.modularity(partition, G)

    dynamic_history.append({
        "r": w,
        "clustering": C,
        "transitivity": TR,
        "assortativity": AS,
        "modularity": M,
        "rewired edges": 0,
        "rewiring_fraction": 0,
        "edge_persistence": 1,
    })

    n = 1  # timestep within stable window
    t = 1  # global timestep
    Ntot = N * N  # for system-spanning flag

    while t <= T and w < f:
        edges = list(G.edges())  # cached for this stable window

        # MAIN DYNAMICS LOOP (stable window)
        while n <= stable_time and t <= T:
            if n % 10 == 0:
                pct = 100 * t / T
                print(f"\r[f={f}, TH={TH}] {pct:6.4f}%  ({t}/{T})", end="", flush=True)

            max_info = max(G.nodes[node]['info'] for node in nodes)
            Delta_range = 1 - max_info

            for node in nodes:
                G.nodes[node]['info'] += random.uniform(0, Delta_range)

            # No avalanche
            if max(G.nodes[node]['info'] for node in nodes) < 1:
                total_info = sum(G.nodes[node]['info'] for node in nodes)
                total_info_transmitted = 0

                history.append({
                    't': t,
                    'avalanche_size': 0,
                    'avalanche_duration': None,
                    'total_info': total_info,
                    'total_info_transmitted': total_info_transmitted,
                    'max_wave_size': 0,
                    'total_topplings': 0,
                    'system_spanning_flag': 0,
                    'branching_mean': None,
                    'branching_std': None,
                    'branching_max': None,
                })

                n += 1
                t += 1
                continue

            # AVALANCHE CASCADE
            avalanche_duration = 0
            total_active_nodes = set()
            currently_active_nodes = []
            previously_active_nodes = []

            max_wave_size = 0
            total_topplings = 0
            cbrs_this_avalanche = []

            while max(G.nodes[node]['info'] for node in nodes) >= 1:
                avalanche_duration += 1
                previously_active_nodes = currently_active_nodes[:]

                currently_active_nodes = [node for node in nodes if G.nodes[node]['info'] >= 1]

                wave_size = len(currently_active_nodes)
                total_topplings += wave_size
                if wave_size > max_wave_size:
                    max_wave_size = wave_size

                if len(previously_active_nodes) > 0:
                    CBR = wave_size / len(previously_active_nodes)
                    branching_history[t].append(CBR)
                    cbrs_this_avalanche.append(CBR)

                total_active_nodes.update(currently_active_nodes)

                for k in currently_active_nodes:
                    I_k = G.nodes[k]['info']
                    NN = list(G.neighbors(k))
                    N_NN = len(NN)
                    if N_NN > 0:
                        share = (A / N_NN) * I_k
                        for nn in NN:
                            G.nodes[nn]['info'] += share
                            G.edges[k, nn]['info_transmitted'] += share
                    G.nodes[k]['info'] = 0

            total_info = sum(G.nodes[node]['info'] for node in nodes)
            total_info_transmitted = sum(G.edges[e]['info_transmitted'] for e in edges)
            avalanche_size = len(total_active_nodes)

            if cbrs_this_avalanche:
                b_mean = sum(cbrs_this_avalanche) / len(cbrs_this_avalanche)
                b_var = sum((x - b_mean) ** 2 for x in cbrs_this_avalanche) / len(cbrs_this_avalanche)
                b_std = b_var ** 0.5
                b_max = max(cbrs_this_avalanche)
            else:
                b_mean = None
                b_std = None
                b_max = None

            if avalanche_size >= 0.9 * Ntot:
                system_spanning_flag = 1
                print("Avalanche spanned system")
            else:
                system_spanning_flag = 0

            history.append({
                't': t,
                'avalanche_size': avalanche_size,
                'avalanche_duration': avalanche_duration,
                'total_info': total_info,
                'total_info_transmitted': total_info_transmitted,
                'max_wave_size': max_wave_size,
                'total_topplings': total_topplings,
                'system_spanning_flag': system_spanning_flag,
                'branching_mean': b_mean,
                'branching_std': b_std,
                'branching_max': b_max,
            })

            n += 1
            t += 1
            if t > T:
                break

        if t > T:
            break

        # Adaptive-network rewiring step.
        unused_edges = [e for e in edges if G.edges[e]['info_transmitted'] < TH]
        n_e = 0

        for (u, v) in unused_edges:
            if random.random() < p:
                n_e += 1
                G.remove_edge(u, v)

                while True:
                    new_v = (random.randint(0, N - 1), random.randint(0, N - 1))
                    if new_v != u and not G.has_edge(u, new_v):
                        G.add_edge(u, new_v)
                        G.edges[u, new_v]['rewiring_binary'] = 1
                        break

        # refresh and reset for next window
        edges = list(G.edges())
        for e in edges:
            G.edges[e]['info_transmitted'] = 0

        m = len(edges)
        rewiring_fraction = n_e / m if m > 0 else 0.0
        persistent_edges = [e for e in edges if G.edges[e]['rewiring_binary'] == 0]
        edge_persistence = (len(persistent_edges) / m) if m > 0 else 0.0

        n = 1
        w += 1

        C = nx.average_clustering(G)
        TR = nx.transitivity(G)
        AS = nx.degree_assortativity_coefficient(G)
        partition = community_louvain.best_partition(G)
        M = community_louvain.modularity(partition, G)

        dynamic_history.append({
            "r": w,
            "clustering": C,
            "transitivity": TR,
            "assortativity": AS,
            "modularity": M,
            "rewired edges": n_e,
            "rewiring_fraction": rewiring_fraction,
            "edge_persistence": edge_persistence,
        })

    parameters = [N, p, T, A, f, TH]
    return history, dynamic_history, branching_history, pos, parameters



# --- summarise dynamic_history into chosen stats ---
def dyn_summary(dh):
    # pull series
    C  = [d["clustering"]      for d in dh]
    TR = [d["transitivity"]    for d in dh]
    AS = [d["assortativity"]   for d in dh]
    M  = [d["modularity"]      for d in dh]
    Er = [d["rewired edges"]   for d in dh]
    Fr = [d["rewiring_fraction"] for d in dh]
    Pi = [d["edge_persistence"]  for d in dh]

    def peak_dev(x):
        x0 = x[0]
        return float(max(abs(v - x0) for v in x))

    def stability(x):
        # mean absolute step change
        return float(np.mean([abs(x[i+1] - x[i]) for i in range(len(x)-1)])) if len(x) > 1 else 0.0

    out = {}
    for name, series in [("clust", C), ("trans", TR), ("assort", AS), ("mod", M)]:
        out[f"{name}_peak_dev"] = peak_dev(series)
        out[f"{name}_stability"] = stability(series)
        out[f"{name}_final"] = float(series[-1])

    out["rewired_edges_mean"] = float(np.mean(Er))
    out["rewired_edges_max"]  = float(np.max(Er))
    out["rewiring_frac_mean"] = float(np.mean(Fr))
    out["rewiring_frac_max"]  = float(np.max(Fr))

    out["edge_persistence_final"] = float(Pi[-1])
    return out


# Sweep
N_list  = [50]
p_list  = [0.2]
T_list  = [100000]
A_list  = [0.5]
f_list  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
TH_list = [0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75]

out_dir = "results/phase_i_d"
os.makedirs(out_dir, exist_ok=True)
summary_path = os.path.join(out_dir, "phase_i_d_dynamic_parameter_sweep.csv")

s_min = 2.5
sys_frac_threshold = 0.9

fields = [
    "N","T","p","A","f","TH","s_min",
    "n_avalanches","n_tail","tau_hat","ks","size_p99","size_max","frac_system_spanning",
    "clust_peak_dev","clust_stability","clust_final",
    "trans_peak_dev","trans_stability","trans_final",
    "assort_peak_dev","assort_stability","assort_final",
    "mod_peak_dev","mod_stability","mod_final",
    "rewired_edges_mean","rewired_edges_max",
    "rewiring_frac_mean","rewiring_frac_max",
    "edge_persistence_final",
]

if not os.path.exists(summary_path):
    with open(summary_path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=fields).writeheader()


# Evaluate each rewiring frequency/threshold combination against the baseline.
param_grid = itertools.product(N_list, p_list, T_list, A_list, f_list, TH_list)

for (N, p, T, A, f, TH) in tqdm(list(param_grid), desc="Running SOC grid"):

    print(f"Running Phase I.D with f={f}, TH={TH}")

    history, dynamic_history, branching_history, pos, parameters = simulate_dynamic_network(N, p, T, A, f, TH)

    sizes = [h["avalanche_size"] for h in history if h.get("avalanche_size", 0) > 0]
    n_aval = len(sizes)

    if n_aval == 0:
        tau_hat, n_tail = np.nan, 0
        ks = np.nan
        size_p99 = np.nan
        size_max = 0
        frac_sys = 0.0
    else:
        size_arr = np.asarray(sizes, dtype=int)
        size_p99 = float(np.percentile(size_arr, 99))
        size_max = int(size_arr.max())
        Ntot = N * N
        frac_sys = float(np.mean(size_arr >= sys_frac_threshold * Ntot))
        tau_hat, n_tail = mle_tau_fixed_xmin(size_arr, s_min)
        ks = ks_distance_fixed_xmin(size_arr, s_min, tau_hat)

    dyn = dyn_summary(dynamic_history)

    row = {
        "N": N, "T": T, "p": p, "A": A, "f": f, "TH": TH, "s_min": s_min,
        "n_avalanches": n_aval, "n_tail": n_tail, "tau_hat": tau_hat, "ks": ks,
        "size_p99": size_p99, "size_max": size_max, "frac_system_spanning": frac_sys,
        **dyn
    }

    with open(summary_path, "a", newline="") as f_out:
        csv.DictWriter(f_out, fieldnames=fields).writerow(row)

    del history, dynamic_history, branching_history

print("Done. Summary saved to:", summary_path)
