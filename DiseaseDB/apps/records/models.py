from django.db import models
from apps.case.models import Patient
from apps.core.models import Record
import uuid
import datetime
import os


def patient_directory_path(instance, filename):
    ext = filename.split('.')[-1].lower()
    new_name = f"{uuid.uuid4().hex[:16]}.{ext if ext in ['jpg', 'png'] else 'jpg'}"
    return os.path.join("compressed", datetime.now().strftime("%Y/%m"), new_name)

class MedicalImage(models.Model, Record):
    # 只存储文件路径，而非文件内容
    compressed_image = models.ImageField(upload_to=patient_directory_path, blank=True)
    patient = models.ForeignKey(Patient, on_delete=models.CASCADE, null=False, related_name='image_set')
    extracted_data = models.JSONField(default=dict)


class TestItem(models.Model):
    RESULT_TYPES = [('quant', '定量'), ('qual', '定性'),]
    item_name = models.CharField(primary_key=True, max_length=100)
    result_type = models.CharField(max_length=10, choices=RESULT_TYPES)
    unit = models.CharField(max_length=20, null=True, blank=True)
    ref_range = models.CharField(max_length=20, null=True, blank=True)

