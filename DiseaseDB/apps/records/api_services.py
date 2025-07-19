import requests
from django.conf import settings


def call_multimodal_api(image_path, json_template):
    url = "https://api.deepseek.com/v1/multimodal"  # 示例API
    headers = {"Authorization": f"Bearer {settings.API_KEY}"}

    # 构建多模态请求
    payload = {
        "prompt": f"提取医疗数据并填充JSON: {json_template}",
        "image": open(image_path, "rb")
    }

    response = requests.post(url, headers=headers, files=payload)
    return response.json()["structured_data"]