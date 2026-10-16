"""Phase II.S: dual-contagion sweep on a static network (HPC workload)."""

try:
    import networkx as nx
except ImportError as e:
    raise SystemExit("networkx is not available on this Python") from e
import numpy as np
from networkx import grid_graph  
import itertools
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **kwargs):
        return x
import os, csv
import random
from scipy.stats import spearmanr

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
def summary_stats(Arr):
    arr = np.asarray([x for x in Arr if x is not None], dtype=float)
    # Remove NaNs if any
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
          return np.nan, np.nan, np.nan
    mean_val = arr.mean()
    p99_val = np.percentile(arr, 99)
    max_val = arr.max()
    return mean_val, p99_val, max_val
def simulate_dual_contagion_static_network(N, p, T, A, B):  # B for the Fomo Contagion Parameter
    # --- Build lattice ---
    G = grid_graph(dim=[N, N])
    pos = {(x, y): (x, y) for x, y in G.nodes()}
    nodes = list(G.nodes())  # small speed win

    # --- Add diagonal edges ---
    for i in range(N):
        for j in range(N):
            if i < N - 1 and j < N - 1:
                G.add_edge((i, j), (i + 1, j + 1))
            if i < N - 1 and j > 0:
                G.add_edge((i, j), (i + 1, j - 1))

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
    total_info = sum(G.nodes[n]['info'] for n in nodes)

    history.append({
        't': 0,
        'avalanche_size': 0,
        'D_avalanche_size': 0,
        'P_avalanche_size': 0,
        'avalanche_duration': None,
        'total_info': total_info,
        'max_D_wave_size': 0,
        'max_P_wave_size': 0,
        'max_wave_size': 0,
        'P_D_ratio': None,
        'P_share_ratio': None,
        'total_topplings': 0,
        'system_spanning_flag': 0,
    })

    n = 1
    Ntot = N * N  # for system-spanning flag

    #                       MAIN DYNAMICS LOOP
    while n < T:
        # --- Compute max info ---
        max_info = max(G.nodes[node]['info'] for node in nodes)
        Delta_range = 1 - max_info

        # --- Update all nodes' info ---
        for node in nodes:
            G.nodes[node]['info'] += random.uniform(0, Delta_range)

        # --- Check for threshold crossing ---
        if max(G.nodes[node]['info'] for node in nodes) < 1:
# Note, we do not update the history if there are no avalanches since we do not care about time between avalanches
            n += 1
            continue

        #                     AVALANCHE CASCADE STARTS
        avalanche_duration = 0

        total_D_active_nodes = set()
        total_P_active_nodes = set()

        # global "already activated (D or P) in this avalanche"
        activated_this_avalanche = set()

        # Per-avalanche bookkeeping
        max_wave_size = 0
        max_D_wave_size = 0
        max_P_wave_size = 0

        total_topplings = 0
        total_D_topplings = 0
        total_P_topplings = 0

        while True:
            currently_active_D_nodes = [node for node in nodes if G.nodes[node]['info'] >= 1]
            if not currently_active_D_nodes:
                break

            avalanche_duration += 1

            active_D_set = set(currently_active_D_nodes)

            # Wave-local (keeps"not both D and P in same wave")
            wave_activated = set(active_D_set)

            # mark D as activated for the whole avalanche
            activated_this_avalanche.update(active_D_set)

            # Wave stats for D
            D_wave_size = len(active_D_set)
            total_D_topplings += D_wave_size
            max_D_wave_size = max(max_D_wave_size, D_wave_size)

            # Collect P candidates (neighbors of D-topplers)
            P_candidates = set()

            total_D_active_nodes.update(active_D_set)

            # topple each active node in D wave
            for k in active_D_set:
                I_k = G.nodes[k]['info']
                NN = list(G.neighbors(k))
                N_NN = len(NN)
                if N_NN > 0:
                    share = (A / N_NN) * I_k
                    for nn in NN:
                        G.nodes[nn]['info'] += share
                    P_candidates.update(NN)
                G.nodes[k]['info'] = 0

            # --- P wave selection: exactly one P-wave attached to this D-wave ---
            currently_active_P_nodes = set()

            for k in P_candidates:
                # never recruit someone twice within the same avalanche
                if k in activated_this_avalanche:
                    continue
                # "not both D and P in same wave"
                if k in wave_activated:
                    continue

                NN = list(G.neighbors(k))
                N_NN = len(NN)
                if N_NN == 0:
                    continue

                active_NN = sum(1 for nb in NN if nb in active_D_set)
                if active_NN == 0:
                    continue

                Prob_Activate = (active_NN / N_NN) ** B
                if random.random() <= Prob_Activate:
                    currently_active_P_nodes.add(k)

            wave_activated.update(currently_active_P_nodes)

            # mark P as activated for whole avalanche
            activated_this_avalanche.update(currently_active_P_nodes)

            # Apply P-wave dynamics
            for k in currently_active_P_nodes:
                I_k = G.nodes[k]['info']
                NN = list(G.neighbors(k))
                N_NN = len(NN)
                if N_NN > 0:
                    share = (A / N_NN) * I_k
                    for nn in NN:
                        G.nodes[nn]['info'] += share

                # partial reset, but if still supercritical, hard-relief it
                if G.nodes[k]['info'] >= 1:
                    G.nodes[k]['info'] = 0
                else:
                    G.nodes[k]['info'] *= 0.5

            total_P_active_nodes.update(currently_active_P_nodes)

            # Wave stats for P
            P_wave_size = len(currently_active_P_nodes)
            total_P_topplings += P_wave_size
            max_P_wave_size = max(max_P_wave_size, P_wave_size)

            wave_size = D_wave_size + P_wave_size
            total_topplings += wave_size
            max_wave_size = max(max_wave_size, wave_size)


        # after avalanche ends
        total_info = sum(G.nodes[node]['info'] for node in G.nodes())
        Dcount = len(total_D_active_nodes)
        Pcount = len(total_P_active_nodes)
        avalanche_size = Dcount + Pcount
        P_D_ratio = (Pcount / Dcount) if Dcount > 0 else None
        P_share_ratio = (Pcount / avalanche_size) if avalanche_size > 0 else None  

        system_spanning_flag = 1 if avalanche_size >= 0.9 * Ntot else 0

        history.append({
            't': n,
            'avalanche_size': avalanche_size,
            'D_avalanche_size': Dcount,
            'P_avalanche_size': Pcount,
            'avalanche_duration': avalanche_duration,
            'total_info': total_info,
            'max_D_wave_size': max_D_wave_size,
            'max_P_wave_size': max_P_wave_size,
            'max_wave_size': max_wave_size,
            'P_D_ratio': P_D_ratio,
            'P_share_ratio': P_share_ratio,
            'total_topplings': total_topplings,
            'system_spanning_flag': system_spanning_flag,
        })
        n += 1

    parameters = [N, p, T, A, B]
    return history
if __name__ == "__main__":
    N = 50
    T = 100000
    p_list = [0.2, 0.5, 0.8] # the slices of the heatmap - 3
    A_list = [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95] # 7
    B_list = [0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3] #11
    out_dir = "results/phase_ii_s"
    os.makedirs(out_dir, exist_ok=True)
    summary_path = os.path.join(out_dir, "phase_ii_s_dual_contagion_sweep.csv")
    s_min = 2.5          
    sys_frac_threshold = 0.9 
    fields = [
        "N","T","p","A","B","s_min",
        "n_avalanches","n_tail","tau_hat","ks",
        "size_p99","size_max","frac_system_spanning", 
        "P_share_mean", "P_share_p99", "P_share_max", "P_ratio_mean", "P_ratio_p99", "P_ratio_max", 
        "frac_P_positive", "spearman_size_Pshare", "spearman_D_P"
    ]
    if not os.path.exists(summary_path):
        with open(summary_path, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=fields).writeheader()
    with open(summary_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)

        total_jobs = len(p_list) * len(A_list) * len(B_list)
        # Sweep topology and both contagion parameters on the static network.
        param_iter = itertools.product(p_list, A_list, B_list)
    
        for (p, A, B) in tqdm(param_iter, total=total_jobs, desc="Running SOC grid"):
            history = simulate_dual_contagion_static_network(N, p, T, A, B)
            avalanches = [h for h in history if h.get("avalanche_size", 0) > 0]
            sizes = [h["avalanche_size"] for h in avalanches]
            P_shares = [h["P_share_ratio"] for h in avalanches]  
            P_ratios = [h["P_D_ratio"] for h in avalanches]    
            n_aval = len(sizes)
            if n_aval == 0:
                tau_hat, n_tail = np.nan, 0
                ks = np.nan
                size_p99 = np.nan
                size_max = 0
                frac_sys = 0.0
                P_share_mean = np.nan 
                P_share_p99 = np.nan 
                P_share_max = np.nan
                P_ratio_mean = np.nan 
                P_ratio_p99 = np.nan
                P_ratio_max = np.nan
                frac_P_positive = np.nan
                spearman_size_Pshare = np.nan
                spearman_D_P = np.nan
            else:
                size_mean, size_p99, size_max = summary_stats(sizes)
                Ntot = N * N
                size_arr = np.asarray(sizes, dtype=float)
                frac_sys = float(np.mean(size_arr >= sys_frac_threshold * Ntot))
                tau_hat, n_tail = mle_tau_fixed_xmin(size_arr, s_min)
                ks = ks_distance_fixed_xmin(size_arr, s_min, tau_hat)
                P_share_mean, P_share_p99, P_share_max = summary_stats(P_shares)
                P_ratio_mean, P_ratio_p99, P_ratio_max = summary_stats(P_ratios)
                S = np.asarray([h["avalanche_size"] for h in avalanches], dtype=float)
                D = np.asarray([h["D_avalanche_size"] for h in avalanches], dtype=float)
                P = np.asarray([h["P_avalanche_size"] for h in avalanches], dtype=float)
                frac_P_positive = float(np.mean(P > 0))
                P_share = np.divide(P, S, out=np.zeros_like(P, dtype=float), where=S > 0)
                spearman_size_Pshare = spearmanr(S, P_share).correlation if len(S) >= 2 else np.nan
                spearman_D_P = spearmanr(D, P).correlation if len(D) >= 2 else np.nan
            row = dict(
                N=N, 
                T=T, 
                p=p, 
                A=A, 
                B=B,
                s_min=s_min, 
                n_avalanches=n_aval, 
                n_tail=n_tail,  
                tau_hat=tau_hat, 
                ks=ks,         
                size_p99=size_p99, 
                size_max=size_max, 
                frac_system_spanning=frac_sys,
                P_share_mean = P_share_mean, 
                P_share_p99 = P_share_p99, 
                P_share_max = P_share_max,
                P_ratio_mean = P_ratio_mean, 
                P_ratio_p99 = P_ratio_p99, 
                P_ratio_max = P_ratio_max,
                frac_P_positive = frac_P_positive,
                spearman_size_Pshare = spearman_size_Pshare,
                spearman_D_P = spearman_D_P,
            )
            with open(summary_path, "a", newline="") as f_out:
                csv.DictWriter(f_out, fieldnames=fields).writerow(row)
            del history
