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
    oligo_concentration = models.PositiveIntegerField(default=0.25e-6)
    monovalent_concentration = models.PositiveIntegerField(default=50e-3)

class SequenceOutput(models.Model):
    sequence = models.OneToOneField(SequenceInput, on_delete=models.CASCADE, primary_key=True)
    output_oligos = models.TextField(max_length=8000)
    primer_dimers = models.TextField(max_length=8000)
    assembly_scheme = models.TextField(max_length=8000, default="Not provided")
    secondary_structure = models.TextField(max_length=8000)
    all_results_url = models.URLField( max_length=300)
