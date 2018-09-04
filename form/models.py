from django.db import models

# Create your models here.


class Contact(models.Model):
    name = models.CharField(max_length=60)
    email = models.CharField(max_length=60)
    title = models.CharField(max_length=60)
    message = models.TextField()

    def __str__(self):
        return self.title