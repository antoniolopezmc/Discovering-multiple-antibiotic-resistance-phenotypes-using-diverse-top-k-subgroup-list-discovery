import pandas as pd
from subgroups.quality_measures import WRAcc
from subgroups.quality_measures import WRAccOptimisticEstimate1
from subgroups.algorithms import VLSD
from datetime import datetime

if __name__ == '__main__':
    # Read the dataset.
    df = pd.read_csv("../mimic-iii-for-experiments.csv")
    print("Number of instances: " + str(len(df)))
    print("Number of attributes: " + str(len(df.columns)))
    target = ("culture_microorganism_name_AND_susceptibility", "ENTEROCOCCUS_SP.-R")
    TP = (df[target[0]] == target[1]).sum()
    FP = (df[target[0]] != target[1]).sum()
    print("INFO: number of instances in which the target is true: " + str(TP))
    print("INFO: number of instances in which the target is false: " + str(FP))
    # Run the VLSD algorithm.
    start_time = datetime.now()
    print("Start time: " + str(start_time))
    model = VLSD(quality_measure = WRAcc(), q_minimum_threshold  = 0.01,
                 optimistic_estimate = WRAccOptimisticEstimate1(), oe_minimum_threshold = 0.01,
                 sort_criterion_in_s1 = VLSD.SORT_CRITERION_NO_ORDER, sort_criterion_in_other_sizes = VLSD.SORT_CRITERION_NO_ORDER,
                 vertical_lists_implementation = VLSD.VERTICAL_LISTS_WITH_BITSETS, write_results_in_file = True, file_path = "subgroups.txt")
    model.fit(df, target)
    print("Number of subgroups mined: " + str(model.selected_subgroups))
    end_time = datetime.now()
    print("End time: " + str(end_time))
