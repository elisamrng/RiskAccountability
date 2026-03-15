import pandas as pd


#Returns the dataframe without repeated activities
def clean_event_log_consecutive(df):
    df = df.copy()
    df = df.sort_values(["case:ID", "time:timestamp"])
    cleaned_parts = []
    affected = []

    for case_id, g in df.groupby("case:ID"):
        seq = g["activity"].tolist()
        keep = [True]
        for i in range(1, len(seq)):
            keep.append(seq[i] != seq[i - 1])
        sub = g.loc[keep].copy()
        if len(sub) < len(g):
            affected.append(case_id)
        cleaned_parts.append(sub)

    df_clean = pd.concat(cleaned_parts, ignore_index=True)
    return df_clean, affected
