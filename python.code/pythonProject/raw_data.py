def rawdata(a1):
    print("原始数据的列名：", a1.columns)
    print("原始数据量：", len(a1))
    t =  a1['遗传模式'].value_counts()
    print ("遗传模式的种类及数量: ",  t)
    # 从a1中选择特定的列并创建新的DataFrame，按照counts列进行排序
    a1 = a1[['gene', 'counts', 'c_change','遗传模式']]
    print("原始数据的前5行：", a1.head())
    #找到基因列表series
    genes_counts = a1['gene'].value_counts()
    print('genes的数量是:' ,len(genes_counts))


    #a2表是整理后的简单表，genes是所有基因的列表
    a2 = a1.copy()
    a2['value'] = a1.apply(lambda row: row['counts'] /33104/2 if row['遗传模式'] == '常染色体' else row['counts'] /33104, axis=1)
    a2 = a2.sort_values(by=['value'], ascending=False)
    print("原始数据包含基因频率value列的前5行：", a2.head())



    # 按照 'gene' 列进行分组
    a2 = a2.groupby('gene')

    # 创建一个字典，用于存储每组的 count 值
    dic = {gene: group['value'].values for gene, group in a2}
    return  dic
