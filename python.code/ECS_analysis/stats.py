import numpy as np
import pandas as pd


class Result:
    def __init__(self, value, interval_list, case, control):
        self.value = value
        self.interval = interval_list
        self.case = case
        self.control = control


def calculate_odds_ratio(input_df, control_label, case_label, gene, label_column):
    # 从input_df提取case，control
    df_control = input_df[input_df[label_column] == control_label]
    df_case = input_df[input_df[label_column] == case_label]

    # 计算四格表中元素abcd
    c = df_control[gene].sum()
    d = df_control.shape[0] - c
    a = df_case[gene].sum()
    b = df_case.shape[0] - a

    # 计算OR值及其95%CI
    odds_ratio = (a * d) / (b * c)
    se_ln = np.sqrt((1/a + 1/b + 1/c + 1/d))
    ci = [np.exp(np.log(odds_ratio)-1.96*se_ln), np.exp(np.log(odds_ratio)+1.96*se_ln)]
    result = Result(odds_ratio, ci, case_label, control_label)
    return result
