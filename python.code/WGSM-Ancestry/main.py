from simu_data import DataGenerator
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    gene = DataGenerator()
    data = gene.generate_geography()
    fig, ax = plt.subplots()

    true_labels = np.repeat(np.arange(0, 5), 1000)
    sns.scatterplot(x=data[:, 0], y=data[:, 1], hue=true_labels)
    plt.show()