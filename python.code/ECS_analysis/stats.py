import numpy as np


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

    '''
    计算四格表中元素x1,x2,n1,n2
              carrier   not_carrier total
    -----------------------
    expose      x1      n1-x1       n1
    non_expose  x2      n2-x2       n2
    '''
    x2 = df_control[gene].sum() + 1
    n2 = df_control.shape[0] + 2
    x1 = df_case[gene].sum() + 1
    n1 = df_case.shape[0] + 2

    # 计算RR值及其95%CI
    rr_value = (x1 * n2) / (n1 * x2)
    se_ln = np.sqrt((1/x1 + 1/x2 - 1/n1 - 1/n2))
    ci = [np.exp(np.log(rr_value)-1.96*se_ln), np.exp(np.log(rr_value)+1.96*se_ln)]
    result = Result(rr_value, ci, case_label, control_label)
    return result
