from django.shortcuts import render
from .forms import ContactForm
from .models import Contact
from django.shortcuts import redirect
from django.urls import reverse
from .forms import ResultRetrievalForm

# Create your views here.

def index(request):
    retrieval_form = ResultRetrievalForm()
    context = {"retrieval_form": retrieval_form}
    return render(request, "form/index.html", context)

def retrieve_results(request):
    sequence_id = ResultRetrievalForm(request.POST)
    if sequence_id.is_valid():
        result_id = request.POST.get("sequence_id")
        redirect_url =  reverse("main:results", args=[result_id])
        return redirect(redirect_url)
    return render(request, "form/contact_form.html")
    # else:
    #     output_oligos = (False, "ERROR56: Invalid Sequence ID", "The provided sequence ID is valid. A valid sequence ID must be 32 digit long alphanumeric sequecne.")
    #     context = {"output_oligos": output_oligos}
    #     return render(request, "main/error.html", context)
        
    

def contact(request):
    contact = ContactForm()
    context = {"contact": contact}
    return render(request, 'form/contact_form.html', context)

def contactOutput(request):
    contact_information = ContactForm(request.POST)
    if contact_information.is_valid():
       name = request.POST.get("name")
       email = request.POST.get("email")
       title = request.POST.get("title")
       message = request.POST.get("message")
       print("Successfully acquired information....")
       contact1 = Contact()
       contact1.name = name
       contact1.email = email
       contact1.title = title
       contact1.message = message
       contact1.save()  
       print("Successfully saved the contact form....") 
       return redirect("forms:contact_message")
    else:
        return redirect("forms:contact_fail")

def contactMessage(request):
    return render(request, "form/contact_success.html")

def contactFail(request):
    return render(request, "form/contact_fail.html")

def viewDocs(request):
    return render(request, "form/documentation.html")

def viewTutorial(request):
    return render(request, "form/tutorial.html")

def viewReferences(request):
    return render(request, "form/references.html")
    
