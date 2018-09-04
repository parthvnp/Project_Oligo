from django.contrib import admin
from .models import Contact

# Register your models here.

class ContactAdmin(admin.ModelAdmin):
    fields = ("name", "title" , "email", "message")

admin.site.register(Contact, ContactAdmin)
