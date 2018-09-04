from django import forms  
from . import choices

class SequenceForm(forms.Form):

    sequence_name = forms.CharField(required=False, widget=forms.TextInput(attrs={"class": "input is-link inputText", "placeholder": "Optional name for the sequence. No special characters allowed.", "name": "sequence_name"}) )

    sequence = forms.CharField(required=True, widget=forms.Textarea(attrs={"class": "textarea is-link inputText", "name": "sequence", "placeholder": "A DNA or RNA sequence containing A, T, U , G and C only. All RNA sequences will be converted to DNA."}))

    max_Length = forms.ChoiceField(required=True, choices=choices.SIZE_CHOICES, initial="50", widget=forms.Select(attrs={"name": "maxLength", "class": "inputText", "value": 50}))

    melt_Tm = forms.DecimalField(required=True, widget=forms.TextInput(attrs={"class": "input is-link inputText", "name": "melt_Tm", "value": 55}), min_value=50, max_value=65)

    oligo_concentration = forms.DecimalField(required=True, widget=forms.TextInput(attrs={ "name": "oligo_concentration" , "class": "input inputText", "value": 0.25}) , min_value=1.0e-12,max_value=10)

    monovalent_concentration = forms.DecimalField(required=True, widget=forms.TextInput(attrs={ "name": "monovalent_concentration",  "class": "input inputText", "value": 5}), min_value = 1.0e-12, max_value=10)




