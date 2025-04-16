
list1=[]

# 遍历列表,转化列表a,为平均值列表 averages
def average(lst):
    while lst:
        # 计算平均值列表
        averages = []
        for i in range(len(lst)):
            ave = sum(lst[:i + 1]) / len(lst[:i + 1])
            averages.append(ave)

        # 找到平均值列表中的最大值及其索引
        max_value = max(averages)
        max_index = len(averages) - 1 - averages[::-1].index(max_value)

        # 将最大值及其索引添加到 list1
        list1.append(max_index)
        list1.append(max_value)

        # 如果最大值的索引是列表的最后一个元素，结束循环
        if len(lst) == 1:
            break

        # 截断列表，继续处理剩余部分
        lst = lst[max_index + 1:]

    return list1

