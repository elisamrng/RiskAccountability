from collections import Counter, defaultdict
from typing import Dict, List
import os
from datetime import datetime

#Extracts the sequence of activities (concept:name) from a PM4Py trace
#e.g. seq = ["Create_Fine", "Send_Fine", ...]
def _sequence_from_rep_trace(trace) -> List[str]:
    result = []
    for ev in trace:
        result.append(ev["concept:name"])
    return result


#Finds from a sequence of activities:
#- repeated activities
#- consecutive repetitions
#- non-consecutive repetitions
def _analyze_sequence(seq: List[str]) -> dict:
    from collections import Counter, defaultdict

    counts = Counter(seq)
    repeats = {}
    
    for a, c in counts.items():
        if c > 1:
            repeats[a] = c
    positions = {}
    for i in range(len(seq)):
        a = seq[i]
        if a not in positions:
            positions[a] = []
        positions[a].append(i)

    consecutive = {}
    i = 0
    while i < len(seq):
        j = i + 1
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        run_len = j - i
        if run_len > 1:
            a = seq[i]
            if a not in consecutive or run_len > consecutive[a]:
                consecutive[a] = run_len
        i = j

    non_consecutive = set()
    for a, pos in positions.items():
        if len(pos) > 1:
            if any(pos[k + 1] - pos[k] > 1 for k in range(len(pos) - 1)):
                non_consecutive.add(a)

    return dict(
        repeats=repeats,
        consecutive=consecutive,
        non_consecutive=sorted(list(non_consecutive)),
        sequence_len=len(seq),
    )


def check_repetitions_from_variants(
    variants: Dict[str, List],
    min_repeats: int = 2,
    show_sequence: bool = True,
):

    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Output"))
    os.makedirs(output_dir, exist_ok=True)
    out_txt = os.path.join(output_dir, "variant_repeats.txt")

    total_variants = len(variants)
    flagged = 0

    with open(out_txt, "w", encoding="utf-8") as f:
        header = (
            f"# Ripetizioni varianti\n"
            f"# Numero varianti analizzate: {total_variants}\n\n"
        )
        f.write(header)

        for idx, (v_id, traces) in enumerate(variants.items(), start=1):

            if not traces:
                continue
            rep_trace = traces[0]
            seq = _sequence_from_rep_trace(rep_trace)
            info = _analyze_sequence(seq)

            if not info["repeats"] and not info["consecutive"]:
                continue
            if not info["repeats"]:
                continue
            if max(info["repeats"].values(), default=0) < min_repeats:
                continue

            flagged += 1
            f.write(f"## Variante {idx}\n")
            f.write(f"Tracce: {len(traces)}\n")
            if show_sequence:
                f.write("Sequenza: " + " -> ".join(seq) + "\n")
            f.write(
                "Attività ripetute: "
                + ", ".join(f"{a} X {c} volte" for a, c in info["repeats"].items())
                + "\n"
            )
            if info["consecutive"]:
                f.write(
                    "  - Massime ripetizioni consecutive: "
                    + ", ".join(f"{a} X {r} volte)" for a, r in info["consecutive"].items())
                    + "\n"
                )
            if info["non_consecutive"]:
                f.write(
                    "  - Ripetizioni non consecutive: "
                    + ", ".join(info["non_consecutive"])
                    + "\n"
                )
            f.write("\n")
        f.write(f"# Varianti con attività ripetute ≥ {min_repeats}): {flagged}\n")
    return flagged
