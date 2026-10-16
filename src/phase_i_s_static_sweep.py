"""Phase I.S: static-network parameter sweep.

Static traditional-network parameter sweep, originally run locally.
Existing CSV rows are preserved on reruns; use an empty output directory for a fresh sweep.
"""

import csv
import itertools
import os
import random
from collections import defaultdict

import networkx as nx
import numpy as np
from networkx import grid_graph
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


def simulate_static_network(N, p, T, A):
    # --- Build lattice ---
    G = grid_graph(dim=[N, N])
    pos = {(x, y): (x, y) for x, y in G.nodes()}
    nodes = list(G.nodes())  # small speed win

    # --- Rewiring ---
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

    # --- History initialisation ---
    history = []
    branching_history = defaultdict(list)  # keep it (per timestep)
    total_info = sum(G.nodes[n]['info'] for n in nodes)

    history.append({
        't': 0,
        'avalanche_size': 0,
        'avalanche_duration': None,
        'total_info': total_info,

        # Record default diagnostics when no avalanche occurs.
        'max_wave_size': 0,
        'total_topplings': 0,
        'system_spanning_flag': 0,
        'branching_mean': None,
        'branching_std': None,
        'branching_max': None,
    })

    n = 1
    Ntot = N * N  # for system-spanning flag

    #                       MAIN DYNAMICS LOOP
    while n < T:
        if n % 10 == 0:   # adjust frequency if you want
            pct = 100 * n / T
            print(f"\r[p={p}, A={A}] {pct:6.4f}%  ({n}/{T})", end="", flush=True)

        # --- Compute max info ---
        max_info = max(G.nodes[node]['info'] for node in nodes)
        Delta_range = 1 - max_info

        # --- Update all nodes' info ---
        for node in nodes:
            G.nodes[node]['info'] += random.uniform(0, Delta_range)

        # --- Check for threshold crossing ---
        if max(G.nodes[node]['info'] for node in nodes) < 1:
            total_info = sum(G.nodes[node]['info'] for node in nodes)
            history.append({
                't': n,
                'avalanche_size': 0,
                'avalanche_duration': None,
                'total_info': total_info,

                # Record per-step diagnostics.
                'max_wave_size': 0,
                'total_topplings': 0,
                'system_spanning_flag': 0,
                'branching_mean': None,
                'branching_std': None,
                'branching_max': None,
            })
            n += 1
            continue

        #                     AVALANCHE CASCADE STARTS
        avalanche_duration = 0
        total_active_nodes = set()
        currently_active_nodes = []
        previously_active_nodes = []

        # Track avalanche-level diagnostics.
        max_wave_size = 0
        total_topplings = 0
        cbrs_this_avalanche = []

        while max(G.nodes[node]['info'] for node in nodes) >= 1:
            avalanche_duration += 1

            previously_active_nodes = currently_active_nodes[:]

            currently_active_nodes = [
                node for node in nodes
                if G.nodes[node]['info'] >= 1
            ]

            # Record wave-level statistics.
            wave_size = len(currently_active_nodes)
            total_topplings += wave_size
            if wave_size > max_wave_size:
                max_wave_size = wave_size

            # branching ratio only after first stage
            if len(previously_active_nodes) > 0:
                CBR = wave_size / len(previously_active_nodes)
                branching_history[n].append(CBR)     
                cbrs_this_avalanche.append(CBR)      # added: per-avalanche summary

            total_active_nodes.update(currently_active_nodes)

            # topple each active node
            for k in currently_active_nodes:
                I_k = G.nodes[k]['info']
                NN = list(G.neighbors(k))
                N_NN = len(NN)
                if N_NN > 0:
                    share = (A / N_NN) * I_k
                    for nn in NN:
                        G.nodes[nn]['info'] += share
                G.nodes[k]['info'] = 0

        # after avalanche ends
        total_info = sum(G.nodes[node]['info'] for node in nodes)
        avalanche_size = len(total_active_nodes)

        # Summarise branching across avalanche waves.
        if cbrs_this_avalanche:
            b_mean = sum(cbrs_this_avalanche) / len(cbrs_this_avalanche)
            # simple std (no numpy)
            b_var = sum((x - b_mean) ** 2 for x in cbrs_this_avalanche) / len(cbrs_this_avalanche)
            b_std = b_var ** 0.5
            b_max = max(cbrs_this_avalanche)
        else:
            b_mean = None
            b_std = None
            b_max = None

        # Identify avalanches reaching the configured system-spanning threshold.
        if avalanche_size >= 0.9 * Ntot:
            system_spanning_flag = 1
            print("Avalanche spanned system")

        else:
            system_spanning_flag = 0

        history.append({
            't': n,
            'avalanche_size': avalanche_size,
            'avalanche_duration': avalanche_duration,
            'total_info': total_info,
            # Record per-step diagnostics.
            'max_wave_size': max_wave_size,
            'total_topplings': total_topplings,
            'system_spanning_flag': system_spanning_flag,
            'branching_mean': b_mean,
            'branching_std': b_std,
            'branching_max': b_max,
        })

        n += 1

    parameters = [N, p, T, A]
    return history, branching_history, pos, G, parameters


N_list = [25]
p_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
T_list = [50000]
A_list = [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.999]
out_dir = "results/phase_i_s"
os.makedirs(out_dir, exist_ok=True)
summary_path = os.path.join(out_dir, "phase_i_s_extended_redistribution_sweep.csv")

s_min = 2.5  # Fixed tail cutoff used in the final Phase I analysis.
sys_frac_threshold = 0.9  # system-spanning if size >= 0.9*N^2

fields = [
    "N","T","p","A","s_min",
    "n_avalanches","n_tail","tau_hat","ks",
    "size_p99","size_max","frac_system_spanning"
]

if not os.path.exists(summary_path):
    with open(summary_path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=fields).writeheader()


# Cartesian sweep across network topology and redistribution strength.
param_grid = list(itertools.product(N_list, p_list, T_list, A_list))
for (N, p, T, A) in tqdm(param_grid, desc="Running SOC grid"):
    print (f"Running Phase I.S (N=25, T=50,000) with p={p}, A={A}")
    history, branching_history, pos, G, parameters = simulate_static_network(N, p, T, A)
     # avalanche sizes only
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

        # system-spanning fraction (derived)
        Ntot = N * N
        frac_sys = float(np.mean(size_arr >= sys_frac_threshold * Ntot))

        # tail fit + KS
        tau_hat, n_tail = mle_tau_fixed_xmin(size_arr, s_min)
        ks = ks_distance_fixed_xmin(size_arr, s_min, tau_hat)

    row = dict(
        N=N, 
        T=T, 
        p=p, 
        A=A, 
        s_min=s_min, #up to here, preset
        n_avalanches=n_aval, #number of avalanches
        n_tail=n_tail,  #potentially determines washout
        tau_hat=tau_hat, #determines color
        ks=ks,          #determines washout
        size_p99=size_p99, #the 99 percentile of avalanches - to see spikes (or outliers)
        size_max=size_max, #largest avalanche - to see spikes
        frac_system_spanning=frac_sys #fraction of system spanning avalanches (more than 90% of grid involved)
    )
    with open(summary_path, "a", newline="") as f:
        csv.DictWriter(f, fieldnames=fields).writerow(row)

    del history, branching_history, G

print("Done. Summary saved to:", summary_path)
