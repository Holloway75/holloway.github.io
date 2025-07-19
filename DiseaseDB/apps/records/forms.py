from django.forms import ModelForm
from .models import MedicalImage, TestItem


class ImageForm(ModelForm):
    class Meta:
        model = MedicalImage
        fields = "__all__"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # 设置字段为不可修改
        self.fields['modified_by'].disabled = True
        self.fields['created_by'].disabled = True