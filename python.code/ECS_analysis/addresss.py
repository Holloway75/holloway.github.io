import os


__all__ = ['Area_counterparts', 'Area_counterparts2', 'Auto_list', 'Xlink_list', 'get_area_from_id', 'Area_sort_list',
           'province_name_simple_to_full']

import pandas as pd

ID_province = {
    "11": "北京",
    "12": "天津",
    '13': '河北',
    '14': '山西',
    '15': '内蒙古',
    '21': '辽宁',
    '22': '吉林',
    '23': '黑龙江',
    '31': '上海',
    '32': '江苏',
    '33': '浙江',
    '34': '安徽',
    '35': '福建',
    '36': '江西',
    '37': '山东',
    '41': '河南',
    '42': '湖北',
    '43': '湖南',
    '44': '广东',
    '45': '广西',
    '46': '海南',
    '50': '重庆',
    '51': '四川',
    '52': '贵州',
    '53': '云南',
    '54': '西藏',
    '61': '陕西',
    '62': '甘肃',
    '63': '青海',
    '64': '宁夏',
    '65': '新疆',
    '75': 'unknown',
    'PA': 'unknown',
    'TZ': 'unknown',
    '83': '台湾',
    '81': '香港',
    '82': '澳门'
}

Autonomous_region = {'新疆':'新疆维吾尔自治区', '广西':'广西壮族自治区', '宁夏':'宁夏回族自治区', '内蒙古':'内蒙古自治区',
                     "西藏":'西藏自治区'}

Province = ['河北', '台湾', '山西', '辽宁', '吉林', '黑龙江', '江苏', '浙江', '安徽', '福建', '江西', '山东', '河南', '湖北',
            '湖南', '广东', '海南', '四川', '贵州', '云南', '陕西', '甘肃', '青海']

Area_counterparts = {
    '广东': ['广东', '海南'],
    '浙沪': ['浙江', '上海'],
    '川渝': ['四川', '重庆'],
    '京津': ['北京', '天津'],
    '云贵': ['云南', '贵州'],
    '闽台': ['福建', '台湾'],
    '西北': ['陕西', '宁夏', '新疆', '甘肃', '青海', '西藏'],
}
Area_counterparts2 = {
    '华南': ['广东', '广西', '海南'],
    '南方': ['浙江', '上海', '江西', '福建', '台湾', '重庆', '湖北', '云南', '贵州', '四川', '湖南', '江苏'],
    '北方': ['北京', '天津','内蒙古', '吉林','陕西', '宁夏', '新疆', '甘肃', '青海', '西藏', '山东', '辽宁', '黑龙江', '河北',
             '安徽', '山西', '河南'],
}
Area_sort_list = ['西北', '黑龙江', '蒙吉', '山西', '河南', '河北', '山东', '辽宁', '安徽', '京津', '江苏', '渝鄂',
       '浙沪', '湖南', '云贵川', '赣闽台', '华南']

os.chdir('D:\我的坚果云\ECS_mid')

Auto_list = open('gene.autosome.list').read().split()

Xlink_list = open('gene.x_link.list').read().split()

Gene_over_200 = ['ABCG5', 'ACADM', 'ACADS', 'ACADSB', 'AGXT', 'ALPL', 'ATP7B', 'CAPN3', 'CEP290', 'CFTR', 'CYP21A2',
                  'CYP27A1', 'ETFDH', 'G6PC1', 'GAA', 'GALC', 'GJB2', 'HBA1/HBA2', 'HBA2', 'HBB', 'MCPH1', 'MMACHC',
                  'MMUT', 'PAH', 'PKHD1', 'POLG', 'PTS', 'SLC22A5', 'SLC25A13', 'SLC26A4', 'SMN1', 'TYR', 'UGT1A1',
                  'USH2A', 'G6PD']


def province_name_simple_to_full(short_name):
    if short_name in Municipality:
        return short_name + '市'
    elif short_name in Autonomous_region.keys():
        return Autonomous_region[short_name]
    elif short_name in Province:
        return short_name + '省'
    elif short_name == 'unknown':
        return short_name
    else:
        raise ValueError


def get_area_from_id(name, fid, data):
    df1 = data[(data.fid == fid) & (data.name == name)]
    if not df1.shape[0]:
        return "unknown"
    elif df1.shape[0] > 1:
        print(df1)
        raise ValueError
    elif pd.isna(data.loc[df1.index[0], 'id']):
        return 'unknown'
    else:
        id_ = str(data.loc[df1.index[0], 'id'])
        return ID_province[id_[0:2]]
