# Generated by Django 2.0.7 on 2018-09-04 04:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('form', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='contact',
            name='title',
            field=models.CharField(default='Not provided', max_length=60),
            preserve_default=False,
        ),
    ]
