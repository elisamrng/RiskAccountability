import os

# The function builds a ranking of the riskiest variants.
# For simplicity, the ranking is produced by reading directly the output
# of the file variants_alignments.txt and extracting the useful data.
# It extracts for each variant the number of traces, the sequence
# and the estimated risk, returning a list of dictionaries
def _parse_variants_alignments(alignments_txt_path: str):
    variants = []
    current = None

    with open(alignments_txt_path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith("## Variante "):
                if current is not None:
                    variants.append(current)


                parts = line.split()
                vidx = None
                if len(parts) >= 3:
                    try:
                        vidx = int(parts[2])
                    except ValueError:
                        vidx = None

                current = {
                    "variant": vidx,
                    "traces": None,
                    "seq": None,
                    "risk": None,
                }
                continue

            if current is None:
                continue

            if line.startswith("Tracce:"):
                try:
                    n_str = line.split(":", 1)[1].strip()
                    current["traces"] = int(n_str)
                except (IndexError, ValueError):
                    pass
                continue

            if line.startswith("Sequenza:"):
                try:
                    seq_str = line.split(":", 1)[1].strip()
                    current["seq"] = seq_str
                except IndexError:
                    pass
                continue

            if line.startswith("Rischio stimato:"):
                try:
                    risk_str = line.split(":", 1)[1].strip()
                    current["risk"] = float(risk_str)
                except (IndexError, ValueError):
                    pass
                continue

        if current is not None:
            variants.append(current)

    variants = [v for v in variants if v.get("risk") is not None]
    return variants


# Generates a report of the ranking of process variants based on the estimated risk.
# The variants are ordered in descending order by risk.
def dump_variants_ranking(
    alignments_txt_path: str,
    out_txt_path: str = "../Output/variants_ranking.txt",
    top_n: int | None = None,
):
    variants = _parse_variants_alignments(alignments_txt_path)

    variants_sorted = sorted(
        variants,
        key=lambda v: (-v["risk"], v["variant"] if v["variant"] is not None else 10**9),
    )

    if top_n is not None and top_n > 0:
        variants_sorted = variants_sorted[:top_n]

    os.makedirs(os.path.dirname(out_txt_path), exist_ok=True)
    with open(out_txt_path, "w", encoding="utf-8") as f:
        if top_n is not None and top_n > 0:
            f.write(f"### Ranking delle prime {len(variants_sorted)} varianti per rischio stimato ###\n")
        f.write("\n")

        for pos, info in enumerate(variants_sorted, start=1):
            vidx = info.get("variant")
            risk = info.get("risk")
            traces = info.get("traces")
            seq = info.get("seq")

            f.write(
                f"{pos}) Variante {vidx} \n"
                f"Tracce: {traces if traces is not None else '?'}\n"
                f"Rischio: {risk:.6f} \n"
            )
            if seq:
                f.write(f"Sequenza: {seq}\n")
            f.write("\n")

    return variants_sorted
