import os
import pm4py
from collections import defaultdict
from utils.risk_utils import (
    load_sigma_map,
    count_deviations_in_alignment,
    combine_union_risk,
    parse_variant_key,
    collapsing_activities,
    sequence_from_trace,
    family_key,
    has_consecutive
)


#Generate a report where for each process variant it reports:
#  - the log-model alignment calculated
#  - the count of model-only and log-only deviations
#  - the calculation of the aggregated severity and the estimated risk
#
#The function can optionally handle families of variants with consecutive repetitions.
#A family of variants groups all variants that present repeated activities
#consecutively: these variants are collapsed into a single clean variant representing the 
# family
def dump_variants_with_alignments_and_risk(
    model_file_pnml,
    variants,
    outfile_path,
    severity_json_path,
    list_avoid=None,
    trace_threshold=10,
    total_traces: int | None = None,
    return_alignments: bool = False,
    variants_raw: dict | None = None,
    target_activity: str | None = "Payment",
    include_family_block: bool = False
):

    sigma_map = load_sigma_map(severity_json_path)

    net, im, fm = pm4py.read_pnml(model_file_pnml, True)

    family_members = {}
    member_sev = {}
    family_sev_mean = {}
    if include_family_block and variants_raw:
        family_members, member_sev, family_sev_mean = compute_family_stats_from_raw(
            variants_raw, net, im, fm, severity_json_path, target_activity=target_activity
        )
    # ----------------------------------------------------
    total_variants = len(variants)
    alignments_by_variant = {}

    os.makedirs(os.path.dirname(outfile_path), exist_ok=True) 
    skipped_outfile = outfile_path.with_name(outfile_path.stem + "_skipped.txt")

    with open(outfile_path, "w", encoding="utf-8") as f, \
         open(skipped_outfile, "w", encoding="utf-8") as fskip:

        header = (
            f"# Elenco allineamenti\n"
            f"# Numero di varianti: {total_variants}\n"
        )
        if total_traces is not None:
            header += f"# Tracce totali (N): {total_traces}\n"
        header += "\n"

        f.write(header)
        fskip.write(header)

        for vidx, v in enumerate(variants, start=1):
            is_family_variant = False
            traces = variants[v]
            if not traces:
                continue

            avoid = False
            if list_avoid:
                for t in traces:
                    case_id = t._get_attributes().get("ID")
                    if case_id in list_avoid:
                        avoid = True
                        break

            rep_trace = traces[0] 
            seq = [ev["concept:name"] for ev in rep_trace]
            align_info = pm4py.algo.conformance.alignments.petri_net.variants.dijkstra_less_memory.apply(
                rep_trace, net, im, fm
            )

            if return_alignments:
                alignments_by_variant[vidx] = align_info

            cost = align_info.get("cost", "n/a")

            cM, cL = count_deviations_in_alignment(align_info["alignment"])
            n_traces = len(traces)
            aggregate_risk = combine_union_risk(cM, cL, sigma_map)

            if total_traces and total_traces > 0:
                p_variant = n_traces / total_traces
            else:
                p_variant = 1.0

            final_risk =  p_variant * aggregate_risk

            target = fskip if avoid else f 
            target.write(f"## Variante {vidx}\n")
            target.write(f"Tracce: {n_traces}\n")
            target.write("Sequenza: " + " -> ".join(seq) + "\n")
            if n_traces <= trace_threshold:
                case_ids = [t._get_attributes().get("ID") for t in traces]
                target.write("Case ID: " + ", ".join(case_ids) + "\n")

            target.write("Allineamento (log | model):\n")
            for (log_move, model_move) in align_info["alignment"]:
                lm = log_move if log_move not in (None, ">>") else ">>"
                mm = model_move if model_move not in (None, ">>") else ">>"
                target.write(f"  ({lm} | {mm})\n")

            target.write(f"Costo algoritmo PM4Py: {cost}\n")
            target.write(f"Severity aggregata σ: {aggregate_risk:.6f}\n")
            # --- verifica preventiva per quale variante è una famiglia ---
            if include_family_block and variants_raw:
                from utils.risk_utils import family_key, parse_variant_key, has_consecutive
                fkey_tmp = family_key(seq)
                for v_key_raw in variants_raw:
                    seq_raw = parse_variant_key(v_key_raw)
                    if has_consecutive(seq_raw, target=target_activity):
                        if tuple(family_key(seq_raw,)) == tuple(
                                fkey_tmp):
                            is_family_variant = True
                            break


            if total_traces and total_traces > 0:
                if not is_family_variant:
                    target.write(f"probabilità λ: {n_traces}/{total_traces} = {p_variant:.6f}\n")


            if include_family_block and variants_raw:

                if '___sev_cache' not in globals():
                    ___sev_cache = {}
                from utils.risk_utils import family_key, parse_variant_key, has_consecutive

                fkey = family_key(seq)
                coll = []
                for v_key_raw, traces_raw in variants_raw.items():
                    seq_raw = parse_variant_key(v_key_raw)
                    if not has_consecutive(seq_raw, target=target_activity):
                        continue
                    if tuple(family_key(seq_raw)) != tuple(
                            fkey):
                        continue
                    m_tuple = tuple(seq_raw)
                    if m_tuple not in ___sev_cache:
                        rep = traces_raw[0] if traces_raw else None #prende traccia rappr.
                        ___sev_cache[m_tuple] = _severity_of_trace(rep, net, im, fm, sigma_map) if rep else 0.0
                    coll.append((m_tuple, len(traces_raw), ___sev_cache[m_tuple]))

                if coll:
                    is_family_variant = True
                    sum_coll_tr = sum(n for (_m, n, _s) in coll)
                    n_raw_clean = n_traces - sum_coll_tr

                    sev_clean = None
                    if n_raw_clean > 0:
                        found = None
                        for v_key_raw, traces_raw in variants_raw.items():
                            if tuple(parse_variant_key(v_key_raw)) == tuple(fkey):
                                found = traces_raw
                                break
                        if found:
                            m_tuple = tuple(fkey)
                            if m_tuple not in ___sev_cache:
                                rep = found[0] if found else None
                                ___sev_cache[m_tuple] = _severity_of_trace(rep, net, im, fm, sigma_map) if rep else 0.0
                            sev_clean = ___sev_cache[m_tuple]
                        else:
                            n_raw_clean = 0

                    if n_traces > 0:
                        weighted_sum = sum(s * n for (_m, n, s) in coll)
                        if n_raw_clean > 0 and (sev_clean is not None):
                            weighted_sum += sev_clean * n_raw_clean
                        fam_sev_weighted = weighted_sum / n_traces
                    else:
                        fam_sev_weighted = 0.0

                    lambda_fam = (n_traces / total_traces) if (total_traces and total_traces > 0) else 1.0

                    target.write("  Sequenza famiglia: " + " -> ".join(fkey) + "\n")
                    target.write("  Varianti collassate (RAW):\n")
                    for (member, tr_count, m_sev) in sorted(coll, key=lambda t: (-t[1], t[0])):
                        target.write(
                            "    - " + " -> ".join(member)
                            + f"  (tracce: {tr_count})  Severity variante: {m_sev:.6f}\n"
                        )
                    if n_raw_clean > 0 and sev_clean is not None:
                        target.write(
                            "    - " + " -> ".join(fkey)
                            + f"  (tracce: {n_raw_clean})  Severity variante: {sev_clean:.6f}\n"
                        )
                    if total_traces and total_traces > 0:
                        target.write(f"  Probabilità famiglia: {n_traces}/{total_traces} = {lambda_fam:.6f}\n")
                    target.write(f"  Severity media famiglia: {fam_sev_weighted:.6f}\n")
                    target.write(
                        f"  Rischio stimato famiglia: {lambda_fam:.6f} × {fam_sev_weighted:.6f} = {(lambda_fam * fam_sev_weighted):.6f}\n")

            if not is_family_variant:
                target.write(f"Rischio stimato: {final_risk:.6f}\n\n")
            else:
                target.write("\n")

    print(f"[OK] report di allineamento {outfile_path}")
    if return_alignments:
        return alignments_by_variant
    return None


# The function generates a textual report containing only the variant families.
# It can be used to better distinguish, among all variants, which ones actually
# contained repetitions and were collapsed into a family.
# This function is optional and can be useful to get a quick overview of what is happening.
def dump_family_map_from_raw_variants(raw_variants: dict, out_txt_path: str, target_activity: str = None):
    fam_to_members = defaultdict(list)
    fam_collapse_set = defaultdict(set)

    for v_key, traces in raw_variants.items():
        seq = parse_variant_key(v_key)
        if not has_consecutive(seq, target=target_activity):
            continue

        target_set = {target_activity} if target_activity is not None else None
        fkey = family_key(seq)
        acts = collapsing_activities(seq)
        if target_activity is not None:
            acts = {a for a in acts if a == target_activity}
        fam_to_members[fkey].append((v_key, len(traces), sorted(acts)))
        fam_collapse_set[fkey].update(acts)

    os.makedirs(os.path.dirname(out_txt_path), exist_ok=True)
    with open(out_txt_path, "w", encoding="utf-8") as f:
        if not fam_to_members:
            f.write("# Nessuna famiglia con ripetizioni consecutive trovata nelle varianti RAW.\n")
        else:
            for fkey in sorted(fam_to_members.keys(), key=lambda k: (len(k), k)):
                fam_str = " -> ".join(fkey)

                if fam_collapse_set[fkey]:
                    f.write(f"## {fam_str} ##\n")
                else:
                    f.write(f"## {fam_str}\n")

                for v_key, n, acts in sorted(fam_to_members[fkey], key=lambda t: (-t[1], t[0])):
                    f.write(f"  - {v_key}  (tracce: {n})\n")
                f.write("\n")

    print(f"[OK] Mappa delle varianti con attività collassate {out_txt_path}")


def _severity_of_trace(trace, net, im, fm, sigma_map):
    import pm4py
    from utils.risk_utils import count_deviations_in_alignment, combine_union_risk
    align = pm4py.algo.conformance.alignments.petri_net.variants.dijkstra_less_memory.apply(
        trace, net, im, fm
    )
    cM, cL = count_deviations_in_alignment(align["alignment"])
    return float(combine_union_risk(cM, cL, sigma_map))


#The function returns:
# - family_members: { fkey(tuple) : [(member(tuple), n_tracce_member), ...] }
# - member_sev:     { member(tuple) : severity_calcolata_su_un_tracciato_rappresentativo }
# - family_sev_mean:{ fkey(tuple) : media_semplice_delle_severity_dei_membri }

def compute_family_stats_from_raw(
    variants_raw: dict,
    net, im, fm,
    severity_json_path: str,
    target_activity: str | None = "Payment",
):

    from collections import defaultdict
    from utils.risk_utils import parse_variant_key, family_key, has_consecutive
    sigma_map = load_sigma_map(severity_json_path)

    family_members = defaultdict(list)
    member_sev = {}

    for v_key, traces in variants_raw.items():
        seq = parse_variant_key(v_key)
        if not has_consecutive(seq, target=target_activity):
            continue
        fkey =family_key(seq)
        family_members[fkey].append((tuple(seq), len(traces)))

    for v_key, traces in variants_raw.items():
        seq = parse_variant_key(v_key)
        if not has_consecutive(seq, target=target_activity):
            continue
        member = tuple(seq)
        if member in member_sev:
            continue  
        if not traces:
            member_sev[member] = 0.0
            continue
        rep_trace = traces[0]  
        sev = _severity_of_trace(rep_trace, net, im, fm, sigma_map)
        member_sev[member] = float(sev)

    family_sev_mean = {}
    for fkey, members in family_members.items():
        vals = [member_sev.get(m, 0.0) for (m, _n) in members]
        family_sev_mean[fkey] = (sum(vals) / len(vals)) if vals else 0.0

    return family_members, member_sev, family_sev_mean

