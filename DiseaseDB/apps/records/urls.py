from django.urls import path
from . import views


urlpatterns = [
    path('/create/', views.TestCreateView.as_view(), name='create_patient'),
    path('/view/', views.PatientListView.as_view(), name='patient_list'),
    path('/detail/<slug:pk>/', views.PatientDetailView.as_view(), name='patient_detail'),
    path('/delete/<slug:pk>/', views.PatientDeleteView.as_view(), name='patient_delete'),
    path('/update/<slug:pk>/', views.PatientUpdateView.as_view(), name='patient_update'),
]