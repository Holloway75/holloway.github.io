import os

__all__ = [
    'area_counterparts',
    'area_counterparts2',
    'Auto_list',
    'Xlink_list',
    'id_province',
    'province_name_simple_to_full',
    'area_counterparts_for_final'
]

id_province = {
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
    '83': '台湾',
    '81': '香港',
    '82': '澳门'
}

Municipality = ['北京', '天津', '上海', '重庆']

Autonomous_region = {'新疆':'新疆维吾尔自治区', '广西':'广西壮族自治区', '宁夏':'宁夏回族自治区', '内蒙古':'内蒙古自治区',
                     "西藏":'西藏自治区'}

Province = ['河北', '台湾', '山西', '辽宁', '吉林', '黑龙江', '江苏', '浙江', '安徽', '福建', '江西', '山东', '河南', '湖北',
            '湖南', '广东', '海南', '四川', '贵州', '云南', '陕西', '甘肃', '青海']

area_counterparts = {
    '桂琼': ['广西', '海南'],
    '浙沪': ['浙江', '上海'],
    '京津': ['北京', '天津'],
    '闽台': ['福建', '台湾'],
    '陕甘宁': ['陕西', '甘肃', '宁夏'],
    '云贵': ['云南', '贵州'],
    '新藏青': ['新疆', '青海', '西藏'],
}

area_counterparts_for_final = {
    '桂琼': ['广西', '海南'],
    '浙沪': ['浙江', '上海'],
    '闽台': ['福建', '台湾'],
    '甘宁': ['甘肃', '宁夏'],
    '新藏青': ['新疆', '青海', '西藏'],
    '广东': ['广东', '香港', '澳门']
}

area_counterparts2 = {
    'Far South': ['广东', '广西', '海南', '香港', '澳门'],
    'South': ['浙江', '上海', '江西', '福建', '台湾', '重庆', '湖北', '云南', '贵州', '四川', '湖南'],
    'North': ['北京', '天津', '内蒙古', '吉林', '陕西', '宁夏', '新疆', '甘肃', '青海', '西藏', '山东', '辽宁', '黑龙江', '河北',
              '江苏', '安徽', '山西', '河南'],
}
os.chdir('D:\我的坚果云\ECS_mid')

with open('gene.autosome.list', 'r') as f:
    Auto_list = f.read().splitlines()

with open('gene.x_link.list', 'r') as f:
    Xlink_list = f.read().splitlines()


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
        print(short_name)
        raise ValueError



