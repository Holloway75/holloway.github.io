import os
from io import BytesIO
from django.core.files.base import ContentFile
from PIL import Image


def process_image(uploaded_file):
    """
    预处理上传图片：验证格式 → 智能压缩 → 安全保存
    参数: uploaded_file (InMemoryUploadedFile) - Django接收的上传文件对象
    返回: 压缩后的ContentFile对象，可直接保存到ImageField
    """
    # 1. 验证图片格式[1,7](@ref)
    ext = os.path.splitext(uploaded_file.name)[1].lower()
    if ext not in ['.jpg', '.jpeg', '.png']:
        raise ValueError("仅支持JPG/PNG格式")

    # 2. 打开图片并转换格式
    img = Image.open(uploaded_file)
    if img.mode != 'RGB':  # 转换PNG透明背景[4](@ref)
        img = img.convert('RGB')

    # 3. 智能压缩（质量+尺寸双重优化）[10,12](@ref)
    output_buffer = BytesIO()
    quality = 85  # 初始质量值
    target_size = 1000 * 1024  # 1MB目标大小

    while quality >= 20 and output_buffer.tell() < target_size:
        output_buffer.seek(0)  # 重置缓冲区
        output_buffer.truncate(0)

        # 质量压缩[11](@ref)
        img.save(output_buffer, format='JPEG', quality=quality, optimize=True)

        # 若仍超限则缩小尺寸（等比缩放0.9倍）[8](@ref)
        if output_buffer.tell() > target_size and quality > 40:
            new_size = (int(img.width * 0.9), int(img.height * 0.9))
            img = img.resize(new_size, Image.LANCZOS)
            quality = max(quality - 5, 20)  # 同时降低质量

    # 4. 生成安全文件名[7](@ref)
    safe_name = f"compressed_{uploaded_file.name.split('.')[0]}.jpg"

    # 5. 返回可直接保存的ContentFile[6](@ref)
    return ContentFile(output_buffer.getvalue(), name=safe_name)