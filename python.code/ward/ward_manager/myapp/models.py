from django.db import models
from django.utils import timezone
from datetime import timedelta


class Diagnosis(models.Model):
    """疾病诊断模板"""
    name = models.CharField(max_length=50)
    icd_code = models.CharField(max_length=10)

    def __str__(self):
        return f"{self.name} ({self.icd_code})"


class LabTest(models.Model):
    """化验检查记录"""
    name = models.CharField(max_length=50)
    order_time = models.DateTimeField(default=timezone.now)
    result_confirmed = models.DateTimeField(null=True, blank=True)
    parameters = models.JSONField()

    @property
    def status(self):
        return 'confirmed' if self.result_confirmed else 'pending'


class MedicalProcedure(models.Model):
    """医疗处置方案"""
    FREQUENCY_CHOICES = [
        ('QD', '每日一次'),
        ('BID', '每日两次'),
        ('TID', '每日三次'),
    ]

    name = models.CharField(max_length=50)
    frequency = models.CharField(max_length=5, choices=FREQUENCY_CHOICES)
    last_performed = models.DateTimeField(null=True, blank=True)
    dosage = models.CharField(max_length=30)


class LabTestProfile(models.Model):
    """疾病关联的化验配置"""
    diagnosis = models.ForeignKey(Diagnosis, on_delete=models.CASCADE, related_name='required_tests')
    lab_type = models.ForeignKey(LabTest, on_delete=models.CASCADE)
    review_interval = models.DurationField()  # 存储timedelta
    critical_threshold = models.CharField(max_length=50)


class Patient(models.Model):
    """住院病人记录"""
    medical_record_no = models.CharField(max_length=20, unique=True)
    diagnosis = models.ForeignKey(Diagnosis, on_delete=models.PROTECT)
    admission_date = models.DateField(default=timezone.now)
    lab_records = models.ManyToManyField(LabTest)
    ongoing_procedures = models.ManyToManyField(MedicalProcedure)

    def get_todo_list(self):
        """生成当日待办事项"""
        todos = []

        # 未确认的检查结果
        pending_labs = self.lab_records.filter(
            order_time__date=timezone.now().date(),
            result_confirmed__isnull=True
        )
        for lab in pending_labs:
            todos.append(f"Review {lab.name} results")

        # 复查提醒逻辑
        for profile in self.diagnosis.required_tests.all():
            latest_lab = self.lab_records.filter(
                lab_type=profile.lab_type
            ).order_by('-order_time').first()

            if latest_lab and latest_lab.result_confirmed:
                due_date = latest_lab.result_confirmed + profile.review_interval
                if timezone.now() > due_date:
                    todos.append(f"Recheck {profile.lab_type.name} (Due {due_date:%Y-%m-%d})")

        # 处置执行跟踪
        for procedure in self.ongoing_procedures.all():
            if not procedure.last_performed or \
                    procedure.last_performed.date() < timezone.now().date():
                todos.append(f"Perform {procedure.get_frequency_display()} {procedure.name}")

        return todos

    def __str__(self):
        return f"Patient {self.medical_record_no}"
