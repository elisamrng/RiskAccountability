import json
from collections import defaultdict

# Loads the JSON file containing severity (sigma) values for activities
# and returns a mapping.
# Missing or invalid values are converted to 0.0.
def load_sigma_map(path):
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f) #transforms the text into an object
    sigma = {} 
    for task, vals in raw.items(): 
        sigma[task] = {
            "sigma_M": float(vals.get("sigma_M", 0.0) or 0.0),
            "sigma_L": float(vals.get("sigma_L", 0.0) or 0.0),
        }
    return sigma


# The function counts the deviations for each activity in an alignment.
#model-only: (log_move in {None, '>>'}) & (model_move not in {None, '>>'})
#log-only:   (model_move in {None, '>>'}) & (log_move not in {None, '>>'})
#Returns:
#counts_M: dictionary activity -> count of model-only deviations
#counts_L: dictionary activity -> count of log-only deviations
def count_deviations_in_alignment(alignment):
    counts_M = defaultdict(int)
    counts_L = defaultdict(int)

    for (log_move, model_move) in alignment: 
        lm_is_void = (log_move in (None, ">>"))
        mm_is_void = (model_move in (None, ">>"))

        if lm_is_void and not mm_is_void:
            counts_M[model_move] += 1
        elif mm_is_void and not lm_is_void:
            counts_L[log_move] += 1
    return counts_M, counts_L


#The function calculates the aggregated risk of a variant by combining
#the model-only and log-only deviations using a union formula
#(1 - product of complements).
#For each activity, the complement is equivalent to (1 − σ)^n,
#where σ is the severity (of log or of model) and n is the number
#of deviations associated with such activity
def combine_union_risk(counts_M, counts_L, sigma_map, multiplicity=1):
    def contrib(sigma, n):
        if sigma <= 0.0 or n <= 0:
            return 0.0
        return 1.0 - (1.0 - sigma) ** n

    factors = []
    for task, n in counts_M.items():
        sigma = sigma_map.get(task, {}).get("sigma_M", 0.0)
        r = contrib(sigma, n)
        factors.append(1.0 - r)

    #Log-only
    for task, n in counts_L.items():
        sigma = sigma_map.get(task, {}).get("sigma_L", 0.0)
        r = contrib(sigma, n)
        factors.append(1.0 - r)

    if not factors:
        return 0.0

    prod = 1.0
    for f in factors:
        prod *= f
    return 1.0 - prod


#The function accepts both strings ('A, B, C') and tuples/lists ('A','B','C').
#It always returns a list of strings.
def parse_variant_key(v_id):
    #Case 1: already a list/tuple -> cast to a list of strings
    if isinstance(v_id, (list, tuple)):
        return [str(x) for x in v_id]

    #Case 2: string -> split on comma (or ';' as fallback)
    if isinstance(v_id, str):
        sep = ',' if ',' in v_id else (';' if ';' in v_id else None)
        if sep:
            return [p.strip() for p in v_id.split(sep) if p.strip()]
        #no separator: single activity
        return [v_id.strip()] if v_id.strip() else []

    #Case 3: any iterable but not string
    try:
        return [str(x) for x in v_id]
    except TypeError:
        return [str(v_id)]


#Returns True if there are consecutive repetitions.
#If target is None: any activity.
#If target is 'Payment' (or another): only that
def has_consecutive(seq, target=None):
    if not seq:
        return False
    prev = seq[0]
    for x in seq[1:]:
        if x == prev and (target is None or x == target):
            return True
        prev = x
    return False

#Returns the activities that appear in consecutive repetitions.
#Example: ['A','A','B','C','C'] -> {'A','C'}
def collapsing_activities(seq):
    res = set()
    if not seq:
        return res
    prev = seq[0]
    for x in seq[1:]:
        if x == prev:
            res.add(x)
        prev = x
    return res


#Given a sequence, collapse consecutive repetitions:
#['A','A','B','B','C'] -> ('A','B','C')
def family_key(seq):
    if not seq:
        return tuple()
    out = [seq[0]]
    for x in seq[1:]:
        if x != out[-1]:
            out.append(x)
    return tuple(out)


#Extract the list of activities from a PM4Py trace
#(uses 'concept:name' or fallback 'activity').
def sequence_from_trace(trace):
    seq = []
    for ev in trace:
        seq.append(ev.get("concept:name", ev.get("activity", "")))
    return seq


#Returns True if the sequence ends with at least two consecutive occurrences
#of the same activity.
def has_trailing_consecutive(seq, target=None):
    if not seq or len(seq) < 2:
        return False
    last = seq[-1]
    if target is not None and last != target:
        return False
    return seq[-2] == last 


#Collapse Family Key that collapses ONLY consecutive repetitions at the end
#Example:
#['A','B','C','C','C'] -> ('A','B','C')
#['A','B','C','C','D'] -> ('A','B','C','C','D') (non collapses in the middle)
#If targets is a set of activities, collapse only if the last one is in targets
def family_key_trailing(seq, targets=None):
    if not seq:
        return tuple()
    last = seq[-1]
    if targets is not None and last not in targets:
        return tuple(seq)
    j = len(seq) - 1
    while j - 1 >= 0 and seq[j - 1] == last:
        j -= 1
    if len(seq) - j < 2:
        return tuple(seq)
    return tuple(seq[:j] + [last])