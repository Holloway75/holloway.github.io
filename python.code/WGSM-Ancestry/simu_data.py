import os

import numpy as np
from scipy.stats import multivariate_normal


class DataGenerator:
    def __init__(self, seed=42):
        np.random.seed(seed)
        self.n_subpops = 5
        self.n_samples = 1000
        self.n_aims = 5000
        self.geo_means = None
        self.geo_covs = None
        self.aim_params = None

    def _generate_geo_distributions(self):
        """生成5个亚群的二维地理分布（示例参数）"""
        # 示例均值和协方差矩阵（可根据实际需求修改）
        self.geo_means = [
            [0, 0], [5, 5], [-5, 5], [5, -5], [-5, -5]
        ]
        self.geo_covs = [
            [[1, 0], [0, 1]],
            [[1, 0.5], [0.5, 1]],
            [[0.5, 0], [0, 0.5]],
            [[2, -0.3], [-0.3, 1]],
            [[1, 0], [0, 2]]
        ]

    def generate_geography(self):
        """生成地理坐标数据"""
        self._generate_geo_distributions()
        coordinates = []
        for mean, cov in zip(self.geo_means, self.geo_covs):
            samples = multivariate_normal.rvs(mean=mean, cov=cov, size=self.n_samples)
            coordinates.append(samples)
        self.geo_data = np.vstack(coordinates)  # shape: (5000, 2)
        return self.geo_data

    def generate_labels(self):
        """生成含噪声的籍贯标签"""
        true_labels = np.repeat(np.arange(self.n_subpops), self.n_samples)

        # 加入5%的标签错误
        error_mask = np.random.rand(len(true_labels)) < 0.05
        error_labels = np.zeros_like(true_labels)

        for i in np.where(error_mask)[0]:
            # 获取当前样本的真实标签
            true_label = true_labels[i]

            # 生成错误标签（确保至少有一个可选标签）
            other_labels = np.delete(np.arange(self.n_subpops), true_label)
            if len(other_labels) == 0:
                # 如果所有标签都被删除（理论上不会发生，但为了鲁棒性）
                error_labels[i] = true_label  # 保持原标签
            else:
                error_labels[i] = np.random.choice(other_labels)

        # 合并正确和错误标签
        noisy_labels = true_labels.copy()
        noisy_labels[error_mask] = error_labels[error_mask]

        self.labels = noisy_labels
        return self.labels

    def _generate_aim_params(self):
        """生成满足条件的AIM位点参数"""
        aim_params = []
        subpop_centers = np.array(self.geo_means)

        for _ in range(self.n_aims):
            while True:
                # 生成参数（可根据地理范围调整正态分布参数）
                a = np.random.normal(0, 0.2) * np.random.choice([-1, 1])
                b = np.random.normal(0, 0.2) * np.random.choice([-1, 1])
                c = np.random.normal(0, 0.5)

                # 计算各亚群中心点的AF
                loci = subpop_centers @ np.array([a, b]) + c
                afs = 1 / (1 + np.exp(-loci))

                # 检查最大AF差异是否在0.15-0.3之间
                max_diff = np.max(afs) - np.min(afs)
                if 0.15 <= max_diff <= 0.3:
                    aim_params.append((a, b, c))
                    break

        self.aim_params = np.array(aim_params)
        return self.aim_params

    def generate_genotypes(self):
        """生成基因型数据"""
        self._generate_aim_params()
        n_total = self.n_subpops * self.n_samples
        genotypes = np.zeros((n_total, self.n_aims), dtype=np.int8)

        # 计算每个样本每个位点的AF
        loci = self.geo_data @ self.aim_params[:, :2].T + self.aim_params[:, 2]
        afs = 1 / (1 + np.exp(-loci))  # shape: (5000, 5000)

        # 生成基因型（基于二项分布）
        for i in range(n_total):
            for j in range(self.n_aims):
                p = afs[i, j]
                genotypes[i, j] = np.random.binomial(2, p)

        self.genotypes = genotypes
        return self.genotypes

    def simulate_sequencing(self, mean_depth=2):
        """模拟低深度测序结果，返回形状为(n_samples, n_aims, 2)的数组
        每个元素包含[总reads数, 突变reads数]，未被覆盖的位点用(-1, -1)表示
        """
        # 初始化三维数组（样本×位点×2）
        seq_data = np.full(
            (self.genotypes.shape[0], self.genotypes.shape[1], 2),
            fill_value=-1,
            dtype=np.int16
        )

        # 生成覆盖掩码（20%覆盖率）
        coverage_mask = np.random.rand(*self.genotypes.shape) < 0.2

        # 生成总reads数（泊松分布）
        total_reads = np.random.poisson(
            lam=mean_depth,
            size=self.genotypes.shape
        ).astype(np.int16)

        # 仅保留被覆盖位点的reads数
        total_reads = np.where(coverage_mask, total_reads, 0)

        # 计算突变等位频率（基因型0→0%，1→50%，2→100%）
        alt_freq = self.genotypes.astype(np.float32) / 2

        # 生成突变reads数（基于二项分布）
        alt_reads = np.random.binomial(
            n=total_reads,
            p=alt_freq
        ).astype(np.int16)

        # 填充数据并处理未覆盖位点
        seq_data[:, :, 0] = total_reads
        seq_data[:, :, 1] = alt_reads
        seq_data[total_reads == 0] = (-1, -1)  # 覆盖但reads数为0视为未覆盖

        self.sequencing_data = seq_data
        return self.sequencing_data


    def save_data(self, filename="simulated_data.npz"):
        """保存所有生成数据（结构化数组需要特殊处理）"""
        np.savez_compressed(
            filename,
            geo_data=self.geo_data,
            labels=self.labels,
            aim_params=self.aim_params,
            genotypes=self.genotypes,
            sequencing_data=self.sequencing_data.view(np.uint8),  # 转换为原始字节
            sequencing_dtype=self.sequencing_data.dtype.descr  # 保存数据类型描述
        )

    def save_data(self, filename="simulated_data.npz"):
        """保存所有生成数据"""
        np.savez_compressed(
            filename,
            geo_data=self.geo_data,
            labels=self.labels,
            aim_params=self.aim_params,
            genotypes=self.genotypes,
            sequencing_data=self.sequencing_data
        )
if __name__ == "__main__":
    os.chdir(r"D:/")
    dg = DataGenerator()
    dg.generate_geography()
    dg.generate_labels()
    dg.generate_genotypes()
    dg.simulate_sequencing()
    dg.save_data()