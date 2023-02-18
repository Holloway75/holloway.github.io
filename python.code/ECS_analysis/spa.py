import pandas as pd
import numpy as np
from scipy.special import comb


def loss_reads(input_df, x1, x2):
    l = 0
    for i in input_df.index:
        a1, a2, b = input_df.loc[i, 'a1'], input_df.loc[i, 'a2'], input_df.loc[i, 'b']
        dp, ad_ref = input_df.loc[i, 'dp'], input_df.loc[i, 'ad_ref']
        fi = freq(a1,a2,b,x1,x2)
        fi_1 = 1 - fi
        if ad_ref == dp:
            c2 = 0.5 ** dp
            l += fi_1 ** 2 + 2 * fi_1 * fi * c2
        elif ad_ref == 0:
            c2 = 0.5 ** dp
            l += 2 * fi * fi_1 * c2 + fi ** 2
        else:
            c2 = comb(dp, ad_ref) * 0.5 ** dp
            l += 2 * fi * fi_1 * c2
    return l


def get_dx_hx_reads(input_df, loc):
    dx, hx = 0, 0
    for i in input_df.index:
        a = np.array([input_df.loc[i , 'a1'], input_df.loc[i, 'a2']]).reshape(2, 1)
        b, dp, ad_ref = input_df.loc[i, 'b'], input_df.loc[i, 'dp'], input_df.loc[i, 'ad_ref']
        fi = freq(a, loc, b)
        fi_1 = 1 - fi
        fi_cross = fi * fi_1
        if ad_ref == dp:
            c2 = 0.5 ** dp
            dx += 2*(-fi_1 + c2 * (1-2*fi)) * fi_cross * a
            hx += -2 * (1 - 4 * fi + 3 * fi2 - c2 * (1 - 6 * fi + 6 * fi2)) * fi_cross * np.dot(a, a.T)
        elif ad_ref == 0:
            c2 = 0.5 ** dp
            dx += 2 * (c2*(1-2*fi)+fi) * fi_cross * a
            hx += -2 * (-c2*(1-6*fi+6*fi2) - (2*fi-3*fi2)) * fi_cross * np.dot(a,a.T)
        else:
            c2 = comb(dp,ad_ref) * 0.5 ** ad_ref * 0.5 ** (dp - ad_ref)
            dx += 2 * (c2*(1-2*fi)) * fi_cross * a
            hx += 2 * c2 * (1-6*fi+6*fi2) * fi_cross * np.dot(a, a.T)
    return [dx, hx]


def get_hx_reads(input_df, loc):
    l = 0
    for i in input_df.index:
        a = np.array([input_df.loc[i , 'a1'], input_df.loc[i, 'a2']]).reshape(2, 1)
        b, dp, ad_ref = input_df.loc[i, 'b'], input_df.loc[i, 'dp'], input_df.loc[i, 'ad_ref']
        fi = freq(a, loc, b)
        fi_1 = 1 - fi
        fi2 = fi ** 2
        if ad_ref == dp:
            c2 = 0.5 ** dp
            l += -2 * (1-4*fi+3*fi2 - c2*(1-6*fi+6*fi2)) * fi * fi_1 * np.dot(a, a.T)
        elif ad_ref == 0:
            c2 = 0.5 ** dp
            l += -2 * (-c2*(1-6*fi+6*fi2) - (2*fi-3*fi2)) * fi * fi_1 * np.dot(a,a.T)
        else:
            c2 = comb(dp,ad_ref) * 0.5 ** ad_ref * 0.5 ** (dp - ad_ref)
            l += 2 * c2 * (1-6*fi+6*fi2) * fi * fi_1 * np.dot(a, a.T)
    return l


