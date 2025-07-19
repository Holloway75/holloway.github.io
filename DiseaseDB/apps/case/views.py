from django.views.generic import CreateView, DeleteView, ListView, DetailView, UpdateView
from django.urls import reverse_lazy
from django.contrib.auth.mixins import LoginRequiredMixin
from guardian.shortcuts import assign_perm, get_objects_for_user
from guardian.mixins import PermissionRequiredMixin
from .models import Patient
from .forms import PatientForm


class PatientCreateView(LoginRequiredMixin, CreateView):
    model = Patient
    form_class = PatientForm
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


class PatientDeleteView(LoginRequiredMixin, PermissionRequiredMixin, DeleteView):
    model = Patient
    template_name = "case/patient_delete_confirm.html"
    success_url = '/patient/view/'
    permission_required = 'case.patient_delete'

    def get_permission_object(self):
        return self.get_object()


class PatientListView(LoginRequiredMixin, ListView):
    model = Patient
    paginate_by = 10
    template_name = "case/patient_list.html"
    context_object_name = 'Patients'

    def get_queryset(self):
        objects = get_objects_for_user(
            user=self.request.user,
            perms=['case.patient_view'],
            klass=Patient,
            any_perm=True,
            with_superuser=True
        )
        return objects


class PatientDetailView(LoginRequiredMixin, PermissionRequiredMixin, DetailView):
    model = Patient
    template_name = 'case/patient_detail.html'
    context_object_name = 'sample'
    permission_required = 'case.patient_view'

    def get_permission_object(self):
        return self.get_object()


class PatientUpdateView(LoginRequiredMixin, PermissionRequiredMixin, UpdateView):
    model = Patient
    form_class = PatientForm
    template_name = 'case/patient.html'
    success_url = reverse_lazy('patient_list')
    permission_required = 'case.patient_change'

    def get_permission_object(self):
        return self.get_object()