from django.forms import ModelForm
from .models import Patient


class PatientForm(ModelForm):
    class Meta:
        model = Patient
        fields = "__all__"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # 设置字段为不可修改
        self.fields['modified_by'].disabled = True
        self.fields['created_by'].disabled = True