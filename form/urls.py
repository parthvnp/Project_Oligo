from django.urls import path, include
from . import views
app_name = "forms"

urlpatterns = [
    path("", views.index, name="input"),
    path("retrive_result/<int: result_id>", views.retrieve_results, name="retrive_result"),
    path("contact", views.contact, name="contact"),
    
]
