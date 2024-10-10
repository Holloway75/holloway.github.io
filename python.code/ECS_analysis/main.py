import os
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["KERAS_BACKEND"] = "jax"
import numpy as np
import pandas as pd
import scipy.stats as stats
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from addresss import *
from ecs_process import AutoGene, XlinkGene, SubPopulation

# def poisson_mixture_analysis(data, gene):
#     columns = ['n_components', 'parameters', 'weights', 'bic']
#     df = pd.DataFrame(columns=columns)
#     for i in np.arange(1,8):
#         mx = PoissonMixture(n_components=i).fit(data)
#         df.loc[i, 'n_components'] = mx.n_components
#         df.loc[i, 'parameters'] = ":".join(mx.params.reshape(-1).astype('str'))
#         df.loc[i, 'weights'] = ":".join(mx.weights.reshape(-1).astype('str'))
#         df.loc[i, 'bic'] = mx.bic(data)
#     if gene == 'HBA1/HBA2':
#         df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\HBA1_HBA2.poisson_mixture.csv', encoding='utf_8_sig')
#     else:
#         df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\%s.poisson_mixture.csv' % gene, encoding='utf_8_sig')


# def poisson(lambda_, k, weight):
#     lambda_, k, weight = float(lambda_), float(k), float(weight)
#     log_prob = -lambda_ + k * np.log(lambda_) - sum([np.log(i) for i in np.arange(1, k+1)])
#     return weight * 15089 * np.exp(log_prob)


# def fac(n):
#     return np.log(np.arange(1, n+1)).sum()


# def binominal(p, k, weight):
#     p, k, weight = float(p), float(k), float(weight)
#     log_prob = fac(1000) - fac(k) - fac(1000-k) + k * np.log(p) + (1000-k) * np.log(1-p)
#     return weight * 15089 * np.exp(log_prob)


# def get_sample_name(dir):
#     os.chdir(dir)
#     fam_list = os.listdir()
#     df = pd.DataFrame(columns=['fid', 'id_1', 'name_1', 'id_2', 'name_2'])
#     for i in fam_list:
#         os.chdir(dir + '\\' + i)
#         sam_list = os.listdir()
#         df.loc[i, 'fid'] = sam_list[0].split("-")[0]
#         df.loc[i, 'id_1'] = sam_list[0].split("-")[1]
#         df.loc[i, 'name_1'] = sam_list[0].split("-")[2]
#         if len(sam_list) == 1:
#             pass
#         elif len(sam_list) == 2:
#             df.loc[i, 'id_2'] = sam_list[1].split("-")[1]
#             df.loc[i, 'name_2'] = sam_list[1].split("-")[2]
#         else:
#             print(i)
#             raise ValueError
#     return df


if __name__ == '__main__':
    os.chdir('D:\我的坚果云\投稿-ecs地区差异')
    df = pd.read_excel('Table 3.xlsx', index_col=0, sheet_name='Sheet2')
    col_list = ['Far South', 'South', 'North']

    # 定义计算地区样本量的数据源和方法
    def get_indiv(_area, _df):
        return _df.loc['individuals_total', _area], _df.loc['individuals_female', _area]

    SubPopulation.data_source = df
    SubPopulation.get_individuals = get_indiv

    # 导入r包
    stats_package = importr('stats')

    # 创建保存结果的df
    gene_list = [i for i in df.index if i in Auto_list + Xlink_list]
    df_result = pd.DataFrame(columns=col_list + ['chi2', 'p', 'OR', 'type'], index=gene_list,
                             data=np.zeros((len(gene_list), 7)))

    for gene_name in gene_list:
        # 创建gene实例
        if gene_name in Auto_list:
            gene = AutoGene(gene_name)
        elif gene_name in Xlink_list:
            gene = XlinkGene(gene_name)
        else:
            raise ValueError
        gene.add_carriers(col_list, df.loc[gene.symbol, col_list])

        # 计算期望频数，用于判断是否需要使用校正或Fisher确切概率法
        carriers = gene.carriers
        col_sum = gene.samples

        df_result.loc[gene.symbol, 'type'] = gene.gene_type
        df_result.loc[gene.symbol, col_list] = carriers

        data = np.array([carriers, col_sum - carriers])
        row_sum = data.sum(axis=1)
        total = col_sum.sum()
        expected = np.outer(row_sum, col_sum) / total

        # 检查期望频数是否都大于5
        min_expected = np.min(expected)

        if min_expected >= 5:
            # 使用普通的卡方检验
            chi2, p, dof, expected = stats.chi2_contingency(data)
            df_result.loc[gene.symbol, ['chi2', 'p']] = [chi2, p]
        else:
            # 使用Fisher确切概率法
            r_data = robjects.r.matrix(robjects.FloatVector(data.flatten()), nrow=2, byrow=True)
            fisher_test = stats_package.fisher_test(r_data)
            p = fisher_test.rx2('p.value')[0]
            df_result.loc[gene.symbol, ['chi2', 'p']] = ['N/A', p]

        # 计算最大OR值，使用拉普拉斯平滑
        ratio = (carriers+1) / (col_sum - carriers + 1 )
        max_OR = np.max(ratio) / np.min(ratio)
        df_result.loc[gene.symbol, 'OR'] = max_OR

    df_result.to_excel('tmp.xlsx', index=True)
