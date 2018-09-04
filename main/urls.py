from django.contrib import admin
from django.urls import path, include
from . import views
app_name = "main"

urlpatterns = [
    path("", views.formInput, name="index"),
    path("results/", views.formOutput, name="output"),
    path("results/view/<uuid:result_id>/", views.viewResults, name="results"),
    path("error/", views.formError, name="error"),
    path("results/download/<uuid:sequence_id>", views.downloadResults, name="download"),
    
]
