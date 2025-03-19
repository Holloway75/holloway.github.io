from django.db import models
from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType
from ward_manager.patients.models import Patient
from django.core.exceptions import ValidationError


class BasicTest(models.Model):
    """基础检验项目"""
    RESULT_TYPES = (
        ('qualitative', '定性'),
        ('quantitative', '定量'),
    )

    test_id = models.CharField("项目ID", max_length=10, unique=True,
                               help_text="格式：B+3位数字，如B001")
    name = models.CharField("项目名称", max_length=100)
    result_type = models.CharField("结果类型", max_length=20, choices=RESULT_TYPES)
    reference = models.CharField("参考值", max_length=200, blank=True)

    class Meta:
        verbose_name = "基础检验项目"

    def __str__(self):
        return f"{self.test_id} - {self.name}"

    def clean(self):
        # 验证ID格式
        if not self.test_id.startswith('B'):
            raise ValidationError("基础项目ID必须以B开头")
        if not self.test_id[1:].isdigit():
            raise ValidationError("ID格式应为B+3位数字")


class Panel(models.Model):
    """检验项目组合"""
    panel_id = models.CharField("组合ID", max_length=10, unique=True,
                                help_text="格式：P+3位数字，如P001")
    name = models.CharField("组合名称", max_length=100)
    tests = models.ManyToManyField(BasicTest, verbose_name="包含项目")

    class Meta:
        verbose_name = "检验组合"

    def __str__(self):
        return f"{self.panel_id} - {self.name}"

    def clean(self):
        # 验证ID格式
        if not self.panel_id.startswith('P'):
            raise ValidationError("组合ID必须以P开头")
        if not self.panel_id[1:].isdigit():
            raise ValidationError("ID格式应为P+3位数字")


class LabTestOrder(models.Model):
    """检验订单"""
    TEST_TYPE_CHOICES = (
        ('basic', '基础项目'),
        ('panel', '项目组合'),
    )
    STATUS_CHOICES = (
        ('pending', '未确认'),
        ('confirmed', '已确认'),
    )

    patient = models.ForeignKey(Patient, verbose_name="患者", on_delete=models.CASCADE)
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.CharField(max_length=10)
    content_object = GenericForeignKey('content_type', 'object_id')
    order_time = models.DateTimeField("开单时间", auto_now_add=True)
    status = models.CharField("状态", max_length=20, choices=STATUS_CHOICES, default='pending')

    class Meta:
        verbose_name = "检验订单"
        indexes = [
            models.Index(fields=['content_type', 'object_id'])
        ]

    def __str__(self):
        return f"{self.patient} - {self.content_object}"


class TestResult(models.Model):
    """检验结果"""
    order = models.ForeignKey(LabTestOrder, verbose_name="检验订单",
                              on_delete=models.CASCADE, related_name='results')
    test = models.ForeignKey(BasicTest, verbose_name="检验项目", on_delete=models.CASCADE)
    quantitative = models.DecimalField("定量结果", max_digits=10, decimal_places=2,
                                       null=True, blank=True)
    qualitative = models.CharField("定性结果", max_length=50, null=True, blank=True)
    created_at = models.DateTimeField("创建时间", auto_now_add=True)

    class Meta:
        verbose_name = "检验结果"
        unique_together = ('order', 'test')  # 防止重复记录

    def clean(self):
        # 验证结果类型匹配
        if self.test.result_type == 'quantitative' and not self.quantitative:
            raise ValidationError("此项目需要输入定量结果")
        if self.test.result_type == 'qualitative' and not self.qualitative:
            raise ValidationError("此项目需要输入定性结果")

    def __str__(self):
        return f"{self.order} - {self.test}"