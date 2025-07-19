from django.shortcuts import render
from .models import MedicalImage
from .services import process_image
from .api_services import call_multimodal_api
from .utils import generate_json_template
from django.views.generic import CreateView, DeleteView, ListView, DetailView, UpdateView
from django.contrib.auth.mixins import LoginRequiredMixin
from .forms import ImageForm


class ItemCreateView(LoginRequiredMixin, CreateView):
    model = MedicalImage
    form_class = ImageForm(initial={'patient': })
    template_name = 'case/patient.html'  # 模板路径
    success_url = '/patient/view/'  # 成功提交后重定向的URL

    def form_valid(self, form):
        _hospital = self.request.user.hospital
        form.instance.owner = _hospital
        form.instance.created_by = self.request.user
        form.instance.modified_by = self.request.user
        response = super().form_valid(form)
        for perm, _ in self.model._meta.permissions:
            full_perm = f"{self.model._meta.app_label}.{perm}"
            assign_perm(full_perm, _hospital, self.object)
        return response



def upload_medical_image(request):
    if request.method == 'POST':
        # 1. 保存原始图片
        image_model = MedicalImage(original_image=request.FILES['image'])
        image_model.savebn()

        # 2. 压缩并保存
        processed_file = process_image(request.FILES['image'])
        image_model.compressed_image.save(
            f"compressed_{request.FILES['image'].name}",
            processed_file
        )

        # 3. 调用多模态API
        json_template = generate_json_template(request.POST['patient_id'])
        extracted_data = call_multimodal_api(
            image_model.compressed_image.path,
            json_template
        )

        # 4. 更新数据库
        image_model.extracted_data = extracted_data
        image_model.save()

        # 5. 填充关联模型
        update_models(extracted_data)  # 根据数据结构更新Patient/TestItem

    return render(request, 'upload.html')