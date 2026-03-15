import os
import pm4py
from collections import defaultdict
from utils.risk_utils import load_sigma_map

#For each activity computes:
# - freq_M: model-only moves per activity
# - freq_L: log-only moves per activity
# - sync: sync moves per activity
#
#The following severity values come from the JSON file:
# - σ_M(A): model severity for activity A
# - σ_L(A): log severity for activity A
#Are used to compute:
# - probability λ_M(A) = freq_M(A) / (freq_M(A) + sync(A))
# - probability λ_L(A) = freq_L(A) / (freq_L(A) + sync(A))
# - risk_M(A) = σ_M(A) * λ_M(A)
# - risk_L(A) = σ_L(A) * λ_L(A)
def dump_activity_severity_from_raw(
    variants_raw: dict,
    net, im, fm,
    severity_json_path: str,
    out_txt_path: str,
):

    sigma_map = load_sigma_map(severity_json_path)

    freq_M = defaultdict(int)
    freq_L = defaultdict(int)
    sync   = defaultdict(int)

    for _v_key, traces in variants_raw.items():
        if not traces:
            continue
        n_traces = len(traces)
        rep_trace = traces[0]

        align = pm4py.algo.conformance.alignments.petri_net.variants.dijkstra_less_memory.apply(
            rep_trace, net, im, fm
        )

        for (log_move, model_move) in align["alignment"]:
            lm_void = (log_move in (None, ">>"))
            mm_void = (model_move in (None, ">>"))

            if lm_void and not mm_void:
                freq_M[model_move] += n_traces
            elif mm_void and not lm_void:
                freq_L[log_move] += n_traces
            else:
                if log_move and model_move and log_move != ">>" and model_move != ">>" and log_move == model_move:
                    sync[log_move] += n_traces

    results = []
    all_acts = set(freq_M.keys()) | set(freq_L.keys()) | set(sync.keys())
    for act in sorted(all_acts):
        m = freq_M.get(act, 0)
        l = freq_L.get(act, 0)
        s = sync.get(act, 0)

        sigma_M = float(sigma_map.get(act, {}).get("sigma_M", 0.0))
        sigma_L = float(sigma_map.get(act, {}).get("sigma_L", 0.0))

        prob_M = (m / (m + s)) if (m + s) > 0 else 0.0
        prob_L = (l / (l + s)) if (l + s) > 0 else 0.0

        risk_M = sigma_M * prob_M
        risk_L = sigma_L * prob_L

        results.append({
            "activity": act,
            "freq_logOnly": l,
            "freq_modelOnly": m,
            "sync": s,
            "sigma_M": sigma_M,
            "sigma_L": sigma_L,
            "prob_M": prob_M,
            "prob_L": prob_L,
            "risk_M": risk_M,
            "risk_L": risk_L,
        })

    results.sort(key=lambda r: max(r["risk_M"], r["risk_L"]), reverse=True)

    os.makedirs(os.path.dirname(out_txt_path), exist_ok=True)
    with open(out_txt_path, "w", encoding="utf-8") as f:
        f.write("# Severity e rischio per attività\n")
        f.write("# σ_M(A), σ_L(A) da severity_map.json\n")
        f.write("# λ_M(A) = freq_M(A) / (freq_M(A) + sync(A))\n")
        f.write("# λ_L(A) = freq_L(A) / (freq_L(A) + sync(A))\n")
        f.write("# rischio_M = σ_M × λ_M ; rischio_L = σ_L × λ_L\n\n")

        for row in results:
            f.write(f"## {row['activity']}\n")
            f.write(f"  freq_logOnly:   {row['freq_logOnly']}\n")
            f.write(f"  freq_modelOnly: {row['freq_modelOnly']}\n")
            f.write(f"  sync:           {row['sync']}\n")
            f.write(f"  σ_M(A):         {row['sigma_M']:.6f}\n")
            f.write(f"  σ_L(A):         {row['sigma_L']:.6f}\n")
            f.write(f"  λ_M(A):         {row['prob_M']:.6f}\n")
            f.write(f"  λ_L(A):         {row['prob_L']:.6f}\n")
            f.write(f"  rischio_M:      {row['risk_M']:.6f}\n")
            f.write(f"  rischio_L:      {row['risk_L']:.6f}\n\n")

    print(f"[OK] Severity e rischi (log/model) per attività {out_txt_path}")
    return results
