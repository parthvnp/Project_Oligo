from django.shortcuts import render

# Create your views here.

def index(request):
    return render(request, "form/index.html")

def retrieve_results(request, result_id):
    return render(request, "form/retrieve_results.html")


def contact(request):
    return render(request, 'form/contact_form.html')