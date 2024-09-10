import copy
import os
import matplotlib.pyplot as plt
import seaborn as sns
import data_prepare
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import plot
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, accuracy_score, roc_auc_score, roc_curve, auc
from sklearn.decomposition import PCA
import statsmodels.api as sm
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["KERAS_BACKEND"] = "jax"
import numpy as np
from addresss import *
import pandas as pd
from pyckmeans import CKmeans
from scipy.stats import chi2


def poisson_mixture_analysis(data, gene):
    columns = ['n_components', 'parameters', 'weights', 'bic']
    df = pd.DataFrame(columns=columns)
    for i in np.arange(1,8):
        mx = PoissonMixture(n_components=i).fit(data)
        df.loc[i, 'n_components'] = mx.n_components
        df.loc[i, 'parameters'] = ":".join(mx.params.reshape(-1).astype('str'))
        df.loc[i, 'weights'] = ":".join(mx.weights.reshape(-1).astype('str'))
        df.loc[i, 'bic'] = mx.bic(data)
    if gene == 'HBA1/HBA2':
        df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\HBA1_HBA2.poisson_mixture.csv', encoding='utf_8_sig')
    else:
        df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\%s.poisson_mixture.csv' % gene, encoding='utf_8_sig')


def poisson(lambda_, k, weight):
    lambda_, k, weight = float(lambda_), float(k), float(weight)
    log_prob = -lambda_ + k * np.log(lambda_) - sum([np.log(i) for i in np.arange(1, k+1)])
    return weight * 15089 * np.exp(log_prob)


def fac(n):
    return np.log(np.arange(1, n+1)).sum()


def binominal(p, k, weight):
    p, k, weight = float(p), float(k), float(weight)
    log_prob = fac(1000) - fac(k) - fac(1000-k) + k * np.log(p) + (1000-k) * np.log(1-p)
    return weight * 15089 * np.exp(log_prob)


def get_sample_name(dir):
    os.chdir(dir)
    fam_list = os.listdir()
    df = pd.DataFrame(columns=['fid', 'id_1', 'name_1', 'id_2', 'name_2'])
    for i in fam_list:
        os.chdir(dir + '\\' + i)
        sam_list = os.listdir()
        df.loc[i, 'fid'] = sam_list[0].split("-")[0]
        df.loc[i, 'id_1'] = sam_list[0].split("-")[1]
        df.loc[i, 'name_1'] = sam_list[0].split("-")[2]
        if len(sam_list) == 1:
            pass
        elif len(sam_list) == 2:
            df.loc[i, 'id_2'] = sam_list[1].split("-")[1]
            df.loc[i, 'name_2'] = sam_list[1].split("-")[2]
        else:
            print(i)
            raise ValueError
    return df


if __name__ == '__main__':
    os.chdir('E:\我的坚果云\胎盘-早产')
    # 创建一些示例数据
    df = pd.read_excel('stats.xlsx', index_col=0, sheet_name='LR')


    # 绘制分组柱状图
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_axes((0.1, 0.25, 0.8, 0.7))

    fontsize=9
    xlabels = [
        'Premature rupture\n of membrane',
        'Perinatal death',
        'Bronchopulmonary\n dysplasia',
        'White matter\n lesion',
        'Intracranial\n hemorrhage',
        'Apgar score\n <8 at 5 min',
        'Meningitis',
        'Retinopathy'
    ]

    # 使用pivot方法重塑数据
    sns.barplot(data=df, x='outcome', y='Log (LRplacenta/LRbasic )', hue='model', alpha=0.9, linewidth=1,
                     edgecolor='w')
    ax.set_xticks(np.arange(0, 8))
    ax.set_xticklabels(xlabels, fontsize=fontsize, rotation=30, ha='center')
    ax.set_yticklabels(np.arange(0, 7), fontsize=fontsize)
    ax.set_ylabel('$Log(LR)$', fontsize=fontsize)
    ax.spines[['right', 'top']].set_visible(False)

    for x_, (text, value) in enumerate(zip(df['p value'], df['Log (LRplacenta/LRbasic )'])):
        c = 'black'
        if text < 0.05:
            c = 'r'
        if x_ == 0:
            text = "p=%.3f" % text
        else:
            text = "%.3f" % text
        x_ = x_ % 8 + (x_ // 8)/2 - 0.45
        ax.text(x_, value+0.1, text, fontsize=fontsize, c=c)
    plt.legend(title='Extraction Methods of Placental Pathological Features', fontsize=fontsize)
    plt.show()








