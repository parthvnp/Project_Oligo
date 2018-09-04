from django import forms

class ResultRetrievalForm(forms.Form):
    sequence_id = forms.UUIDField(required=True, widget=forms.TextInput(attrs={"class": "input inputText is-fullwidth", "placeholder": "Sequence ID to retrieve results", "name": "sequenceID", "size": "36"}))

class ContactForm(forms.Form):
    name = forms.CharField(required=True, widget=forms.TextInput(attrs={"class": "input is-link", "name": "name", "placeholder": ""}))
    email = forms.EmailField(required=True, widget=forms.EmailInput(attrs={"class": "input is-link", "name": "email",}))
    title = forms.CharField(required=True, widget=forms.TextInput(attrs={"class": "input is-link", "name": "title", "placeholder": ""}))
    message = forms.CharField(required=True, widget=forms.Textarea(attrs={"class": "textarea is-link", "name": "message",  }))