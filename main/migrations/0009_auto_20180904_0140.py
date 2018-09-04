# Generated by Django 2.0.7 on 2018-09-04 01:40

import datetime
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0008_auto_20180904_0133'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='sequenceoutput',
            name='id',
        ),
        migrations.AlterField(
            model_name='sequenceinput',
            name='sequence_name',
            field=models.CharField(default=datetime.datetime(2018, 9, 4, 1, 40, 34, 947731), max_length=60),
        ),
        migrations.AlterField(
            model_name='sequenceoutput',
            name='sequence',
            field=models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, serialize=False, to='main.SequenceInput'),
        ),
    ]
