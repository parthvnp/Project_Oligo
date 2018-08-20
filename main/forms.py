from django import forms  
from . import choices

class SequenceForm(forms.Form):

    sequence_name = forms.CharField(required=False, widget=forms.TextInput(attrs={"class": "input is-link", "placeholder": "optional", "name": "sequence_name"}) )

    sequence = forms.CharField(required=True, widget=forms.Textarea(attrs={"class": "textarea is-link", "name": "sequence", "placeholder": "A DNA or RNA sequence containing A, T, U , G and C only. All RNA sequences will be converted to DNA."}))

    max_Length = forms.ChoiceField(required=True, choices=choices.SIZE_CHOICES, widget=forms.Select(attrs={"name": "maxLength", "default": 50}))

    melt_Tm = forms.DecimalField(required=True, widget=forms.TextInput(attrs={"class": "input is-link", "name": "melt_Tm", "value": 55}), min_value=50, max_value=65)

    algorithm = forms.ChoiceField(choices=choices.ALGORITHM_CHOICES, widget=forms.Select(attrs={"name": "algorithm", "default": "Zuker"}))



