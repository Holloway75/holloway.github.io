from django.db.models.signals import pre_delete
from django.dispatch import receiver
from django.conf import settings
import shutil
import os
from apps.case.models import Patient

@receiver(pre_delete, sender=Patient)
def delete_patient_media(sender, instance, **kwargs):
    """删除患者时同步删除其媒体目录"""
    media_path = os.path.join(settings.MEDIA_ROOT, 'patients', str(instance.patient_id))
    if os.path.exists(media_path):
        shutil.rmtree(media_path)  # 递归删除目录及内容