def generate_json_template(patient_id):
    # 动态关联模型字段
    patient_fields = [f.name for f in Patient._meta.fields if f.name not in ['id']]
    test_item_fields = [f.name for f in TestItem._meta.fields if f.name not in ['id', 'image']]

    return {
        "patient_data": {field: "" for field in patient_fields},
        "test_items": [{field: "" for field in test_item_fields}]
    }