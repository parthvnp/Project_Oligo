from django.db import models
from . import choices
import uuid
import datetime
# Create your models here.


class SequenceInput(models.Model):
    id = models.UUIDField(primary_key=True, editable=False, default=uuid.uuid4)
    sequence_name = models.CharField(max_length=60, default=datetime.datetime.now())
    sequence = models.TextField()
    max_Length = models.PositiveIntegerField(choices=choices.SIZE_CHOICES, default=50)
    melt_Tm = models.PositiveIntegerField(default=55)
    algorithm = models.CharField(choices=choices.ALGORITHM_CHOICES, default="Zuker", max_length=20)

class SequenceOutput(models.Model):
    sequence = models.ForeignKey(SequenceInput, on_delete=models.CASCADE)
    output_oligo = models.TextField(max_length=80)
    oligo_melt_tm = models.PositiveIntegerField()
    secondary_structure = models.TextField(max_length=80)
    secondary_structure_url = models.URLField( max_length=300)
