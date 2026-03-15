import os

# Sorts the activities based on a score computed as a weighted
# combination of the model-only risk and the log-only risk.
# The score is defined as:
# score = delta_model · risk_M + delta_log · risk_L
# Returns the list of activities sorted by decreasing risk.
def rank_linear(metrics, delta: float = 1.0, delta_log: float = 1.0, top_n: int | None = None):
    out = []
    for r in metrics:
        score = delta * r.get('risk_M', 0.0) + delta_log * r.get('risk_L', 0.0)
        e = dict(r)
        e['score'] = score
        out.append(e)
    out.sort(key=lambda x: (x['score'], x.get('risk_M',0.0), x.get('risk_L',0.0)), reverse=True)
    if top_n:
        out = out[:top_n]
    return out

# Sorts without considering a score (a function originally planned but no longer used).
# Sorts the activities lexicographically, using a primary and a secondary
# ordering criterion (for example risk_M and risk_L).
def rank_lexicographic(metrics, primary='risk_M', secondary='risk_L', top_n: int | None = None):
    out = sorted(metrics, key=lambda x: (x.get(primary,0.0), x.get(secondary,0.0)), reverse=True)
    if top_n:
        out = out[:top_n]
    return out


#Computes the activity ranking using rank_linear
#stores the result in a text file
def dump_activity_ranking_linear(metrics, out_txt_path: str, delta: float = 1.0, delta_log: float = 1.0, top_n: int | None = None):
    os.makedirs(os.path.dirname(out_txt_path) or '.', exist_ok=True)
    ranked = rank_linear(metrics, delta=delta, delta_log=delta_log, top_n=top_n)
    with open(out_txt_path, 'w', encoding='utf-8') as f:
        f.write(f"### Ranking delle prime {top_n} attività per rischio stimato ###\n")

        for i, r in enumerate(ranked, start=1):
            f.write(f"{i}) {r['activity']}\n")
            f.write(f"   score:        {r['score']:.6f}\n")
            f.write(f"   rischio_M:    {r.get('risk_M', 0.0):.6f}\n")
            f.write(f"   rischio_L:    {r.get('risk_L', 0.0):.6f}\n")
            f.write(f"   σ_M:          {r.get('sigma_M', 0.0):.6f}\n")
            f.write(f"   σ_L:          {r.get('sigma_L', 0.0):.6f}\n")
            f.write(f"   λ_M:          {r.get('prob_M', 0.0):.6f}\n")
            f.write(f"   λ_L:          {r.get('prob_L', 0.0):.6f}\n")
            f.write(f"   freq_M:       {int(r.get('freq_modelOnly', 0))}\n")
            f.write(f"   freq_L:       {int(r.get('freq_logOnly', 0))}\n")
            f.write(f"   sync:         {int(r.get('sync', 0))}\n\n")
    return ranked

def dump_activity_ranking_lexi(metrics, out_txt_path: str, primary='risk_M', secondary='risk_L', top_n: int | None = None):
    os.makedirs(os.path.dirname(out_txt_path) or '.', exist_ok=True)
    ranked = rank_lexicographic(metrics, primary=primary, secondary=secondary, top_n=top_n)
    with open(out_txt_path, 'w', encoding='utf-8') as f:
        f.write('# Ranking attività\n')
        f.write(f'# Ordine: {primary} (desc), poi {secondary} (desc). top_n={top_n if top_n else "all"}\n\n')
        for i, r in enumerate(ranked, start=1):
            f.write(f"{i}) {r['activity']}\n")
            f.write(f"   rischio_M: {r.get('risk_M',0.0):.6f}   rischio_L: {r.get('risk_L',0.0):.6f}\n")
            f.write(f"   σ_M: {r.get('sigma_M',0.0):.6f}  σ_L: {r.get('sigma_L',0.0):.6f}   λ_M: {r.get('prob_M',0.0)::.6f}  λ_L: {r.get('prob_L',0.0):.6f}\n")
            f.write(f"   freq_M: {int(r.get('freq_modelOnly',0))}   freq_L: {int(r.get('freq_logOnly',0))}   sync: {int(r.get('sync',0))}\n\n")
    return ranked
