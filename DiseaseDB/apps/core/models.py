from django.db import models
from django.contrib.auth import get_user_model
from django.contrib.auth.models import Group


User = get_user_model()


class Record(models.Model):
    """
    抽象基类，用于记录数据创建者、所有者和最后修改者，
    以及创建时间和最后修改时间
    """
    created_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, blank=True,
                                   related_name='%(class)s_created_by', verbose_name='创建者')

    owner = models.ForeignKey(Group, on_delete=models.SET_NULL, null=True, blank=True,
                              related_name='%(class)s_owned_by', verbose_name='所有者')

    modified_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, blank=True,
                                    related_name='%(class)s_modified_by', verbose_name='最后修改者')

    created_at = models.DateTimeField(auto_now_add=True, verbose_name='创建时间')

    modified_at = models.DateTimeField(auto_now=True, verbose_name='最后修改时间')

    class Meta:
        abstract = True

    @property
    def get_owner(self):
        return f"{self.owner}"
