from django.contrib.auth.models import AbstractUser, Group
from django.db import models


class CustomUser(AbstractUser):
    # 添加额外字段
    name = models.CharField(max_length=30, null=True, blank=True, verbose_name='姓名')
    phone = models.CharField(max_length=20, null=True, blank=True, verbose_name="手机号")
    hospital = models.ForeignKey(
        Group,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        verbose_name="医院",
        related_name="subordinate_users"
    )

    # 修改默认字段属性（可选）
    email = models.EmailField(unique=True, verbose_name="电子邮箱")

    def get_name(self):
        return self.name

    class Meta:
        verbose_name = '用户'
        verbose_name_plural = '用户列表'
        db_table = 'custom_users'  # 自定义表名

    def __str__(self):
        return f"{self.name}"