# Risk Analysis in Business Processes

This repository contains the work carried out during the internship of:

**Emanuele Machetti**
Bachelor's Degree in Computer Science
Department of Computer Science
University of Turin

Under the supervision of 
**Matteo Baldoni, Cristina Baroglio, Elisa Marengo, and Roberto Micalizio**

------------------------------------------------------------------------

## Repository Structure

``` text
Dataset/        # Contains event log files (.XES)
PetriNet/       # Contains process models (.PNML)
Scripts/        # Main scripts
 ├─ main.py
 ├─ ...
 └─ utils/      # Utility scripts
Config/         # Configuration files for assigning severity weights to activities (severity_map.json)
Output/         # Contains generated outputs
```

------------------------------------------------------------------------

## Required Packages

Running `main.py` requires the packages listed in `requirements.txt`.

Install them with:

    python -m pip install -r requirements.txt

------------------------------------------------------------------------

## Configuration Guide

The `main.py` script is the main entry point of the project, from which
the entire execution flow and the produced output can be configured.

Before proceeding, make sure you have:

-   an event log file in **.XES** format
-   the corresponding process model in **.PNML** format
-   a description of the activity severities in a **.JSON** format

Log files must be placed inside the **Dataset/** folder, while process
models must be placed inside the **PetriNet/** folder.
severity_map.json must be placed insidse **Config/** folder.

------------------------------------------------------------------------

## Quick Execution

Move to the `Scripts/` folder and run:

    python main.py

------------------------------------------------------------------------

## Generated Outputs

The outputs are text reports stored in the `Output/` folder.

``` text
Output/
 ├─ activity_ranking.txt
 ├─ activity_severity.txt
 ├─ variants_ranking.txt
 ├─ variant_repeats.txt
 ├─ variants_alignments.txt
 └─ variants_alignments_skipped.txt
```


### activity_ranking.txt

Ranking of activities with the highest severity score.

### activity_severity.txt

Lists the severity values associated with each activity.

### variants_ranking.txt

Ranking of the most risky variants.


### variant_repeats.txt

Reports the number of repetitions of each repeated activity within each
variant.

### variants_alignments.txt

Displays the alignment between log traces and the process model together
with the associated severity values.

### variants_alignments_skipped.txt

Contains only the variants that were excluded from the main analysis
using `list_avoid`.
