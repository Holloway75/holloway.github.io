from django.db import models
from django.contrib.auth import get_user_model
from apps.core.models import Record

User = get_user_model()

class Patient(Record, models.Model):
    """
    患者模型，继承自Record抽象基类
    主键是patient_id，所有属性均为选填
    """
    patient_id = models.CharField(primary_key=True, max_length=20, verbose_name='患者ID')
    patient_name = models.CharField(max_length=30, null=True, blank=True, verbose_name='患者姓名')
    GENDER_CHOICES = [
        ('M', '男'),
        ('F', '女'),
        ('O', '其他'),
        ('U', '未知'),
    ]
    gender = models.CharField(max_length=1, choices=GENDER_CHOICES, null=True, blank=True, verbose_name='性别')
    birth_date = models.DateField(null=True, blank=True, verbose_name='出生日期')
    primary_diagnosis = models.CharField(max_length=255, null=True, blank=True, verbose_name='主要诊断')

    class Meta:
        verbose_name = '患者'
        verbose_name_plural = '患者'
        permissions = [
            ('patient_view', 'Can view owned resource'),
            ('patient_change', 'Can change owned resource'),
            ('patient_delete', 'Can delete owned resource'),
        ]

    def __str__(self):
        return f"{self.patient_name} {self.patient_id}"
