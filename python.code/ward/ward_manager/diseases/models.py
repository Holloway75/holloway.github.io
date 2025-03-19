from django.db import models
from ward_manager.patients.models import Patient
from ward_manager.labtests.models import BasicTest, Panel, LabTestOrder
from django.utils import timezone
from datetime import timedelta

class Disease(models.Model):
    """疾病表"""
    name = models.CharField("疾病名称", max_length=100, unique=True)
    description = models.TextField("描述", blank=True)
    related_basic_tests = models.ManyToManyField(
        BasicTest,
        verbose_name="关联基础化验项目",
        blank=True
    )
    related_panels = models.ManyToManyField(
        Panel,
        verbose_name="关联化验组合",
        blank=True
    )
    recheck_interval = models.PositiveIntegerField(
        "复查间隔天数",
        default=7,
        help_text="重要化验项目需要复查的时间间隔（天）"
    )
    important_tests = models.ManyToManyField(
        BasicTest,
        verbose_name="重要化验项目",
        related_name="important_diseases",
        blank=True
    )

    def __str__(self):
        return self.name

class PatientDisease(models.Model):
    """患者疾病关联表"""
    patient = models.ForeignKey(Patient, verbose_name="患者", on_delete=models.CASCADE)
    disease = models.ForeignKey(Disease, verbose_name="疾病", on_delete=models.CASCADE)
    diagnosed_at = models.DateTimeField("确诊时间", auto_now_add=True)
    attending_doctor = models.ForeignKey(
        'users.User',
        verbose_name="主治医师",
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )

    def create_initial_orders(self):
        """创建初始化验订单"""
        # 创建关联的基础项目订单
        for test in self.disease.related_basic_tests.all():
            LabTestOrder.objects.create(
                patient=self.patient,
                content_object=test,
                test_type='basic',
                status='pending'
            )
        # 创建关联的panel订单
        for panel in self.disease.related_panels.all():
            LabTestOrder.objects.create(
                patient=self.patient,
                content_object=panel,
                test_type='panel',
                status='pending'
            )

class TestReminder(models.Model):
    """化验复查提醒表"""
    RESPONSE_CHOICES = (
        ('pending', '待处理'),
        ('postpone', '明日提醒'),
        ('done', '已复查'),
        ('cancel', '不复查')
    )

    patient_disease = models.ForeignKey(
        PatientDisease,
        verbose_name="患者疾病",
        on_delete=models.CASCADE
    )
    test = models.ForeignKey(
        BasicTest,
        verbose_name="化验项目",
        on_delete=models.CASCADE
    )
    due_date = models.DateField("提醒日期")
    last_test_date = models.DateField("上次检验日期", null=True, blank=True)
    status = models.CharField(
        "状态",
        max_length=20,
        choices=RESPONSE_CHOICES,
        default='pending'
    )

    def update_after_test(self, test_date):
        """更新检验后的提醒状态"""
        self.last_test_date = test_date
        self.due_date = test_date + timedelta(days=self.patient_disease.disease.recheck_interval)
        self.status = 'pending'
        self.save()

    def handle_response(self, response):
        """处理医生响应"""
        if response == 'postpone':
            self.due_date += timedelta(days=1)
        elif response == 'cancel':
            self.status = 'cancel'
        elif response == 'done':
            self.status = 'done'
            # 创建结果确认提醒
            ResultConfirmationReminder.objects.create(
                test_order=self.get_last_test_order(),
                due_date=timezone.now().date()
            )
        self.save()

class ResultConfirmationReminder(models.Model):
    """结果确认提醒表"""
    CONFIRM_CHOICES = (
        ('pending', '待处理'),
        ('confirmed', '已确认不录入'),
        ('entered', '已录入结果'),
        ('postpone', '明日提醒')
    )

    test_order = models.ForeignKey(
        LabTestOrder,
        verbose_name="化验订单",
        on_delete=models.CASCADE
    )
    due_date = models.DateField("提醒日期")
    status = models.CharField(
        "确认状态",
        max_length=20,
        choices=CONFIRM_CHOICES,
        default='pending'
    )

    def handle_confirmation(self, response):
        """处理确认响应"""
        if response == 'postpone':
            self.due_date += timedelta(days=1)
        else:
            self.status = response
            if response == 'entered':
                self.test_order.status = 'confirmed'
                self.test_order.save()
        self.save()