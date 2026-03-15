# run with python and not with python3
"""
Process Mining Analysis Tool for Event Log Conformance and Risk Assessment
This module performs comprehensive process mining analysis on event logs against Petri net models.
It includes:
Event Log Processing:
    - Reads event logs from CSV or XES formats
    - Extracts and manages process variants 
Petri Net Construction:
    - Creates simple linear Petri nets from activity sequences
    - Maintains proper start/end markings for conformance checking
Conformance Analysis:
    - Computes optimal alignments between log traces and process model
    - Calculates fitness, precision, and other conformance metrics
    - Identifies deviations: model-only skips (unexecuted model activities)
    - Identifies deviations: log-only skips (unmodeled log activities)
Risk and Severity Assessment:
    - Maps activities to severity levels from JSON configuration
    - Ranks activities and variants by severity and conformance risk
    - Computes aggregate risk metrics per activity and variant
Main Workflow:
    - Input: Petri net model (PNML) + Event log (XES/CSV) + Severity map (JSON)
    - Processing: Clean, align, analyze conformance and risk
    - Output: Multiple text reports covering alignments, rankings
"""
from pathlib import Path

import warnings

warnings.filterwarnings("ignore", message="The argument 'infer_datetime_format' is deprecated")
warnings.filterwarnings("ignore", message="the EventLog class has been deprecated")
warnings.filterwarnings("ignore", category=UserWarning)

import pm4py
import pandas as pd
from utils import input_file
from utils.risk_dump import dump_variants_with_alignments_and_risk
from utils.check_repetitions import check_repetitions_from_variants
from utils.risk_dump import dump_family_map_from_raw_variants
from utils.risk_dump_single_activity import dump_activity_severity_from_raw
from utils.consecutive_clean import clean_event_log_consecutive
from utils.activity_ranking import (
    dump_activity_ranking_linear)
from utils.variants_ranking import dump_variants_ranking

def get_variants_xes(log_file_csv):
    dataframe = pd.read_csv(log_file_csv, sep=';')
    dataframe = dataframe.rename(columns={'ID': 'case:ID'})
    dataframe['activity'] = dataframe['activity'].apply(lambda x: x.replace(' ','_'))
    df = pm4py.format_dataframe(dataframe, case_id='case:ID', activity_key='activity', timestamp_key='time:timestamp')
    event_log = pm4py.convert_to_event_log(df)
    variants = pm4py.statistics.variants.log.get.get_variants(event_log)
    return variants


# Loads an event log and extracts the process variants.
# If USE_CLEAN_ALIGNMENT = True, it applies a cleaning step that removes repeated activities.
# Calls check_repetitions_from_variants to identify variants containing repetitions.
# Returns:
#   - event_log: the possibly cleaned log
#   - variants: variants after the cleaning step
#   - variants_raw: original unmodified variants
def get_variants(log_file, csv=True, single_case = False, case_id = -1):
    if csv: # reads from a csv
        dataframe = pd.read_csv(log_file, sep=';')
        dataframe = dataframe.rename(columns={'ID': 'case:ID'})
        
    else: # reads from xes
        dataframe = pm4py.read_xes(log_file)
        dataframe = dataframe.rename(columns={'case:concept:name': 'case:ID'})
        dataframe = dataframe.rename(columns={'concept:name': 'activity'})

    dataframe['activity'] = dataframe['activity'].apply(lambda x: x.replace(' ','_'))
    df = pm4py.format_dataframe(dataframe, case_id='case:ID', activity_key='activity', timestamp_key='time:timestamp')
    
    if (single_case):
        df = dataframe.loc[dataframe['case:ID'] == case_id]

    event_log_raw = pm4py.convert_to_event_log(df)
    variants_raw = pm4py.statistics.variants.log.get.get_variants(event_log_raw)
    event_log = pm4py.convert_to_event_log(df)
    num_cases = len(event_log)

    if USE_CLEAN_ALIGNMENT:
        print("Cleaning consecutive repetitions...")
        df_tmp = pm4py.convert_to_dataframe(event_log_raw)
        df_clean, affected = clean_event_log_consecutive(df_tmp)
        event_log = pm4py.convert_to_event_log(df_clean)
    else:
        print("Standard Alignment...")
        event_log = event_log_raw

    variants = pm4py.statistics.variants.log.get.get_variants(event_log)
    check_repetitions_from_variants(variants)
    return event_log, variants, variants_raw


# Given a model and a set of variants, returns a dictionary where they keys
# are the activities and the values are the total number of skips for the 
# activities of the optimal alignments of every log trace
def find_activity_skips(model_file_pnml, variants, list_avoid):
    net, initial_marking, final_marking = pm4py.read_pnml(model_file_pnml, True)

    activities = {}    # sum model-only + log-only
    models_skip = {}   # model-only
    logs_skip = {}     # log-only

    i = 0
    for v in variants:
        i += 1
        avoid = any(trace._get_attributes().get("ID") in list_avoid for trace in variants[v])
        if avoid:
            continue

        number_of_traces = len(variants[v])
        alignment = pm4py.algo.conformance.alignments.petri_net.variants.dijkstra_less_memory.apply(
            variants[v][0], net, initial_marking, final_marking
        )

        for (log_move, model_move) in alignment['alignment']:
            # MODEL-ONLY: log is ">>", model is the activity
            if (log_move in (None, ">>")) and (model_move not in (None, ">>")):
                models_skip[model_move] = models_skip.get(model_move, 0) + number_of_traces
                activities[model_move] = activities.get(model_move, 0) + number_of_traces

            # LOG-ONLY: log is the activity, model is ">>"
            elif (model_move in (None, ">>")) and (log_move not in (None, ">>")):
                logs_skip[log_move] = logs_skip.get(log_move, 0) + number_of_traces
                activities[log_move] = activities.get(log_move, 0) + number_of_traces

    return i, activities, models_skip, logs_skip

# This method computes the model path corresponding to the alignment of each trace.
# Returns a dictionary where the key is the model path and the value is a dictionary containing
# the alignments of the traces to the path and the total number of traces aligning to it.
def find_model_path(model_file_pnml, variants, list_avoid):

    net, initial_marking, final_marking = pm4py.read_pnml(model_file_pnml, True) # model net
    model_paths = {}

    print("\nNumber of variants: \n")
    print(f"{len(variants)}\n")

    for v in variants:
        case_list = []
        avoid = False

        for case in variants[v]:
            case_id = case._get_attributes()["ID"]
            case_list.append(case_id)
            if case_id in list_avoid:
                avoid = True
        
        if not avoid:
            number_of_traces = len(variants[v])
            alignment = pm4py.algo.conformance.alignments.petri_net.variants.dijkstra_less_memory.apply(variants[v][0], net, initial_marking, final_marking)
            alignemnt2 = pm4py.algo.conformance.alignments.petri_net.algorithm.apply_trace(variants[v][0], net, initial_marking, final_marking)
            model_path = [model for (log,model) in alignment['alignment'] if None!=model and model != ">>"]
            model_path = tuple(model_path)
            if model_path in model_paths:
                model_paths[model_path]["num_traces"] = model_paths[model_path]["num_traces"] + number_of_traces
                model_paths[model_path]["alignment"].append(alignment['alignment'])
            else:
                model_paths[model_path] = {"num_traces" : number_of_traces, "alignment" : [alignment['alignment']]}
    return model_paths
        

def print_solutions(s):
    for k in s.keys():
        print("Model Trace: \n", k, "\n")
        print("Number of traces: ", s[k]['num_traces'], "\n")
        print("Alignments with Variants: ")
        for a in s[k]['alignment']:
            print(a)
        print("\n\n\n")


if __name__ == '__main__':
    
    BASE_DIR = Path(__file__).resolve().parent.parent
    input_file_pnml_to_clean = str(BASE_DIR / "PetriNet" / "FineManagementModel.pnml")

    #--------------- Severity File ----------------
    severity_path_file = str(BASE_DIR / "Config" / "severity_map.json")
    #---------------------------------------------

    #--------------- Input log File ---------------
    dataset_file_xes = str(BASE_DIR / "Dataset" / "Road_Traffic_Fine_Management_Process.xes")
    #---------------------------------------------

    #--------------- Configuration Parameters ---------------
    USE_CLEAN_ALIGNMENT = False                 
    TOP_K_VARIANTS_FOR_PATTERNS = 10 
    # 
    #---------------------------------------------
    
    model_file_pnml = str(BASE_DIR / "PetriNet" / "model.pnml") 
    dataset_file_csv = str(BASE_DIR / "Dataset" / "Test-short.csv")
    family_txt = str(BASE_DIR / "Output" / "family_severity.txt")
    family_map_txt = str(BASE_DIR / "Output" / "families_map.txt")

    # The input file needs to be cleaned from graphic data
    # Present if the net is produced with GreatSPN
    input_file.clean(input_file_pnml_to_clean, model_file_pnml)

    # the import of the log can be done from csv or from xes file
    is_input_csv = False     
    case_id = "A10001"
    # Some cases can be avoided if needed
    list_avoid = []
    
    net, initial_marking, final_marking = pm4py.read_pnml(model_file_pnml, True)

    fitness = None
    precision = None
    evaluation_values = None
    num_cases = 0

    if is_input_csv:
        event_log, variants, variants_raw = get_variants(dataset_file_csv, is_input_csv)
        fitness = pm4py.fitness_alignments(event_log, net, initial_marking, final_marking)
        precision = pm4py.algo.evaluation.precision.algorithm.apply(event_log, net, initial_marking, final_marking)
        evaluation_values = pm4py.algo.evaluation.algorithm.apply(event_log, net, initial_marking, final_marking)
        num_cases = len(event_log)
    else:
        event_log, variants, variants_raw = get_variants(dataset_file_xes, is_input_csv)
        fitness = pm4py.fitness_alignments(event_log, net, initial_marking, final_marking)
        precision = pm4py.algo.evaluation.precision.algorithm.apply(event_log, net, initial_marking, final_marking)
        evaluation_values = pm4py.algo.evaluation.algorithm.apply(event_log, net, initial_marking, final_marking)
        num_cases = len(event_log)
        print("\n\nNumber of cases: ", num_cases)
        print("\n\nNumber of variants: ", len(variants))

    
    i, activities, models_skip, logs_skip = find_activity_skips(model_file_pnml, variants, list_avoid)

    print("\nAll deviations (model-only + log-only):")
    for a, c in sorted(activities.items(), key=lambda kv: kv[1], reverse=True):
        print(f"  {a}: {c}")

    print("\nModel-only deviations:")
    for a, c in sorted(models_skip.items(), key=lambda kv: kv[1], reverse=True):
        print(f"  {a}: {c}")

    print("\nLog-only deviations:")
    for a, c in sorted(logs_skip.items(), key=lambda kv: kv[1], reverse=True):
        print(f"  {a}: {c}")


    print("\n\nFitness event log: ", fitness)
    print("\n\nPrecision event log: ", precision)
    print("\n\nEvaluation event log dictionary with measures: ", evaluation_values, "\n")
    print("\n\nNumber of Cases: ", num_cases, "\n")

    
    outfile = BASE_DIR / "Output" / "variants_alignments.txt"

    alignments_by_variant = dump_variants_with_alignments_and_risk(
        model_file_pnml, variants,
        outfile_path=outfile,
        severity_json_path=severity_path_file,
        list_avoid=list_avoid,
        trace_threshold=10,
        total_traces=num_cases,
        return_alignments = True,
        variants_raw=variants_raw,
        target_activity=None,
        include_family_block=USE_CLEAN_ALIGNMENT
    )

    if USE_CLEAN_ALIGNMENT:
        # Families map
        dump_family_map_from_raw_variants(
            raw_variants=variants_raw,
            out_txt_path=family_map_txt,
            target_activity=None,
        )

    # --- severity single activities ---
    activity_sev_txt = BASE_DIR / "Output" / "activity_severity.txt"

    metrics = dump_activity_severity_from_raw(
        variants_raw=variants_raw,
        net=net, im=initial_marking, fm=final_marking,
        severity_json_path=severity_path_file,
        out_txt_path=activity_sev_txt,
    )

    # --- activity ranking ---
    activity_ranking = dump_activity_ranking_linear(
        metrics,
        out_txt_path = BASE_DIR / "Output" / "activity_ranking.txt",
        delta=1, #peso rischio_M
        delta_log=1, #peso rischio_L
        top_n=TOP_K_VARIANTS_FOR_PATTERNS
    )
    
    out_txt_path = BASE_DIR / "Output" / "variants_ranking.txt",
    # --- variant ranking
    variants_ranking = dump_variants_ranking(
        alignments_txt_path=outfile,
        out_txt_path = str(BASE_DIR / "Output" / "variants_ranking.txt"),
        top_n=TOP_K_VARIANTS_FOR_PATTERNS
    )