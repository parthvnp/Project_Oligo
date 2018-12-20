from django.urls import path, include
from . import views
app_name = "forms"

urlpatterns = [
    path("", views.index, name="index"),
    path("retrieve_results/", views.retrieve_results, name="retrieve_results"),
    path("contact", views.contact, name="contact"),
    path("contact/output", views.contactOutput, name="contact_output"),
    path("contact/message", views.contactMessage, name="contact_message"),
    path("contact/fail", views.contactFail, name="contact_fail"),
    path("documentation/", views.viewDocs, name="docs"),
    path ("tutorial/", views.viewTutorial, name="tutorial"),
    path("references/", views.viewReferences, name="references"),
    
]
