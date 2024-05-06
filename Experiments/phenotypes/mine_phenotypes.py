import pandas as pd
from subgroups.algorithms import DSLM
from subgroups.data_structures import SubgroupList
from subgroups.core import Subgroup
from subgroups.quality_measures.quality_measure import QualityMeasure
from subgroups.quality_measures.wracc import WRAcc
from bitarray import bitarray
from typing import Union
from datetime import datetime
from re import compile
import sys

def mine_diverse_topk_subgroup_lists(dataset : pd.DataFrame,
                                     subgroups_input_file_path : str,
                                     target_attribute_name : str,
                                     target_attribute_value : Union[str, int, float],
                                     maximum_number_of_subgroup_lists : int,
                                     maximum_number_of_subgroups_per_subgroup_list : int,
                                     beta_parameter : float, 
                                     maximum_positive_overlap : float,
                                     maximum_negative_overlap : float,
                                     diverse_topk_subgroup_lists_file_path) -> None:
    # Run the DSLM algorithm.
    model = DSLM(input_file_path = subgroups_input_file_path,
                 max_sl = maximum_number_of_subgroup_lists,
                 sl_max_size = maximum_number_of_subgroups_per_subgroup_list,
                 beta = beta_parameter,
                 maximum_positive_overlap = maximum_positive_overlap,
                 maximum_negative_overlap = maximum_negative_overlap,
                 output_file_path = diverse_topk_subgroup_lists_file_path
                 )
    model.fit(dataset, (target_attribute_name, target_attribute_value))

def load_diverse_topk_subgroup_lists(dataset : pd.DataFrame,
                                     TP : int,
                                     FP : int,
                                     target_attribute_name : str,
                                     target_attribute_value : str,
                                     diverse_topk_subgroup_lists_file_path : str,
                                     quality_measure : QualityMeasure) -> (list, float, float):
    sl_index = 0
    sl_quality_dict = dict()
    sl_coverage = pd.Series([False]*len(dataset))
    list_of_subgroup_lists = []
    dataset_mask = (dataset[target_attribute_name] == target_attribute_value)
    number_of_dataset_instances = len(dataset)
    input_file = open(diverse_topk_subgroup_lists_file_path, "r")
    regex_object_sl_header = compile("^## Subgroup list (?P<n_subgroups>\(.*\)) ##$")
    regex_object_subgroup = compile("^s(?P<number>[0-9]+): Description: (?P<description>.+), Target: (?P<target>.+)$")
    for line in input_file:
        match_object_1 = regex_object_sl_header.fullmatch(line.rstrip("\n"))
        match_object_2 = regex_object_subgroup.fullmatch(line.rstrip("\n"))
        if match_object_1:
            list_of_subgroup_lists.append( SubgroupList(bitarray(dataset_mask.tolist(), endian = "big"), bitarray((~dataset_mask).tolist(), endian = "big"), number_of_dataset_instances) )
            sl_index = sl_index + 1
            sl_quality_dict[sl_index] = 0
        elif match_object_2:
            subgroup = Subgroup.generate_from_str("Description: " + match_object_2.group("description") + ", Target: " + match_object_2.group("target"))
            tp_Series, fp_Series, _  = subgroup.filter(dataset)
            list_of_subgroup_lists[-1].add_subgroup(subgroup, bitarray(tp_Series.tolist(), endian = "big"), bitarray(fp_Series.tolist(), endian = "big"))
            subgroup_qm = quality_measure.compute({QualityMeasure.TRUE_POSITIVES : tp_Series.sum(), QualityMeasure.FALSE_POSITIVES : fp_Series.sum(), QualityMeasure.TRUE_POPULATION : TP, QualityMeasure.FALSE_POPULATION : FP})
            sl_quality_dict[sl_index] = sl_quality_dict[sl_index] + subgroup_qm
            sl_coverage = sl_coverage | tp_Series | fp_Series
    input_file.close()
    return list_of_subgroup_lists, sum(sl_quality_dict.values())/len(sl_quality_dict.values()), sl_coverage.sum()/len(dataset)

#grid = {
#    "beta_parameter" : [0.0, 0.15, 0.3, 0.4, 0.5, 0.65, 0.75, 0.9, 1.0],
#    "maximum_positive_overlap" : [0.00, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.60],
#    "maximum_negative_overlap" : [0.00, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.60]
#}

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("ERROR: invalid number of parameters.")
        print(str(sys.argv))
        exit(1)
    # Read the program parameters.
    b = float(sys.argv[1])
    pos_o = float(sys.argv[2])
    neg_o = float(sys.argv[3])
    # Create and open the log file.
    log_file = open("log/log_" + str(b) + "_" + str(pos_o) + "_" + str(neg_o) + ".txt", "w")
    # Read the dataset.
    df = pd.read_csv("../mimic-iii-for-experiments.csv")
    print("Number of instances: " + str(len(df)) + "\n", end="")
    log_file.write("Number of instances: " + str(len(df)) + "\n")
    print("Number of attributes: " + str(len(df.columns)) + "\n", end="")
    log_file.write("Number of attributes: " + str(len(df.columns)) + "\n")
    target = ("culture_microorganism_name_AND_susceptibility", "ENTEROCOCCUS_SP.-R")
    TP = (df[target[0]] == target[1]).sum()
    FP = (df[target[0]] != target[1]).sum()
    print("Number of instances in which the target is true: " + str(TP) + "\n", end="")
    log_file.write("Number of instances in which the target is true: " + str(TP) + "\n")
    print("Number of instances in which the target is false: " + str(FP) + "\n", end="")
    log_file.write("Number of instances in which the target is false: " + str(FP) + "\n")
    # Executing the DSLM algorithm.
    print("#########################################################\n", end="")
    log_file.write("#########################################################\n")
    print("#########################################################\n", end="")
    log_file.write("#########################################################\n")
    print("beta_parameter = " + str(b) + "\n", end="")
    log_file.write("beta_parameter = " + str(b) + "\n")
    print("maximum_positive_overlap = " + str(pos_o) + "\n", end="")
    log_file.write("maximum_positive_overlap = " + str(pos_o) + "\n")
    print("maximum_negative_overlap = " + str(neg_o) + "\n", end="")
    log_file.write("maximum_negative_overlap = " + str(neg_o) + "\n")
    output_file_name = "output/" + str(b) + "_" + str(pos_o) + "_" + str(neg_o) + ".txt"
    start_time = datetime.now()
    mine_diverse_topk_subgroup_lists(dataset = df,
                        subgroups_input_file_path  = "input_subgroups.txt",
                        target_attribute_name = target[0],
                        target_attribute_value = target[1],
                        maximum_number_of_subgroup_lists = 3,
                        maximum_number_of_subgroups_per_subgroup_list = 10,
                        beta_parameter = b, 
                        maximum_positive_overlap = pos_o,
                        maximum_negative_overlap = neg_o,
                        diverse_topk_subgroup_lists_file_path = output_file_name)
    end_time = datetime.now()
    list_of_subgroup_lists, quality, coverage = load_diverse_topk_subgroup_lists(df, TP, FP, target[0], target[1], output_file_name, WRAcc())
    print("---\n", end="")
    log_file.write("---\n")
    for sl in list_of_subgroup_lists:
        print(str(sl) + "\n", end="")
        log_file.write(str(sl) + "\n")
    print("---\n", end="")
    log_file.write("---\n")
    print("quality = " + str(quality) + "\n", end="")
    log_file.write("quality = " + str(quality) + "\n")
    print("coverage = " + str(coverage) + "\n", end="")
    log_file.write("coverage = " + str(coverage) + "\n")
    print("---\n", end="")
    log_file.write("---\n")
    print("Start time: " + str(start_time) + "\n", end="")
    log_file.write("Start time: " + str(start_time) + "\n")
    print("End time: " + str(end_time) + "\n", end="")
    log_file.write("End time: " + str(end_time) + "\n")
    # Close the log file.
    log_file.close()
