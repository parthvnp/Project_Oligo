from django.urls import path
from . import views
app_name = "dimercheck"
urlpatterns = [
path("", views.index, name="index"),

]