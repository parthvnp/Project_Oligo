from django.shortcuts import render, redirect
from . forms import SequenceForm
from .models import SequenceInput
from oligodesign import OligoDesign
from django.urls import reverse

def index(request):
    redirect_url = reverse("form", "main.urls")
    return redirect(redirect_url)

def formInput(request):
    form = SequenceForm()
    context = {"form": form}
    return render(request, "main/form_input.html", context)

def formOutput(request):
    form = SequenceForm(request.POST)
    if form.is_valid():
        sequence_name = request.POST.get("sequence_name")
        sequence = request.POST.get("sequence")
        max_Length = request.POST.get("max_Length")
        melt_Tm = request.POST.get("melt_Tm")
        algorithm = request.POST.get("algorithm")
        s1 = SequenceInput()
        s1.sequence_name = sequence_name
        s1.sequence = sequence
        s1.max_Length = max_Length
        s1.melt_Tm = melt_Tm
        s1.algorithm = algorithm
        s1.save()
        redirect_url = reverse("results", "main.urls", args=(s1.pk,))
        return redirect(redirect_url)
    return render(request, "404.html")


def formError(request):
   return render(request, "404.html")


def viewResult(request, result_id):
    s1 = SequenceInput.objects.get(pk=result_id)
    sequence = s1.sequence
    output = OligoDesign.EquiTmOligo(sequence=sequence)
    output_oligos = output.main()
    output_oligo_profile = output.oligos_representation()
    primer_dimer = output.primer_dimer_check()
    context = {"output_oligos": output_oligos, "output_oligos_profile": output_oligo_profile, "result_id": result_id, "primer_dimer": primer_dimer}
    return render(request, "main/results.html", context)