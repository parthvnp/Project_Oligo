from django.contrib import admin
from django.urls import path, include
from . import views
app_name = "main"

urlpatterns = [
    path("", views.formInput, name="index"),
    path("results/", views.formOutput, name="output"),
    path("results/<uuid:result_id>/", views.viewResult, name="results"),
    path("error/", views.formError, name="error"),
    
]
