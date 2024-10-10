import warnings
import numpy as np


class SubPopulation:
    data_source = None
    get_individuals = None          # 传入函数，作为类方法

    def __init__(self, area):
        self.area = area
        self._individuals_total, self._individuals_female = SubPopulation.get_individuals(area, self.data_source)

    @property
    def total(self):
        return self._individuals_total

    @property
    def female(self):
        return self._individuals_female


class Gene:
    def __init__(self, gene_symbol: str):
        self.symbol = gene_symbol
        self.samples = None
        self.area_carriers_dict = {}      # 创建地区-携带者数字典
        self.gene_type = None

    def add_carriers(self, areas, carriers):
        for key, value in zip(areas, carriers):
            if key in self.area_carriers_dict.keys():
                warnings.warn(f'"{key}" was changed in area_carrier_dict')
            self.area_carriers_dict[key] = value

    @property
    def carriers(self):
        return np.array([self.area_carriers_dict[i] for i in self.area_carriers_dict.keys()])


class AutoGene(Gene):
    def __init__(self, gene_symbol):
        super().__init__(gene_symbol)
        self.gene_type = 'Auto'

    def add_carriers(self, areas, carriers):
        super().add_carriers(areas, carriers)
        self.samples = np.array([SubPopulation(area).total for area in self.area_carriers_dict.keys()])


class XlinkGene(Gene):
    def __init__(self, gene_symbol):
        super().__init__(gene_symbol)
        self.gene_type = 'Xlink'

    def add_carriers(self, areas, carriers):
        super().add_carriers(areas, carriers)
        self.samples = np.array([SubPopulation(area).female for area in self.area_carriers_dict.keys()])
