from django.shortcuts import render, redirect
from django.http import HttpResponse
from . forms import SequenceForm
from .models import SequenceInput, SequenceOutput
from oligodesign import OligoDesign
from django.urls import reverse
from django.core.exceptions import ObjectDoesNotExist

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
        oligo_concentration = request.POST.get("oligo_concentration")
        monovalent_concentration = request.POST.get("monovalent_concentration")
        s1 = SequenceInput()
        s1.sequence_name = sequence_name
        s1.sequence = sequence
        s1.max_Length = max_Length
        s1.melt_Tm = melt_Tm
        oligo_concentration = float(oligo_concentration) * 10**-6
        s1.oligo_concentration = oligo_concentration
        monovalent_concentration = float(monovalent_concentration) * 10**-3
        s1.monovalent_concentration = monovalent_concentration
        s1.save()
        redirect_url = reverse("main:results", args=(s1.pk,))
        return redirect(redirect_url)
    return render(request, "404.html")


def formError(request):
   return render(request, "404.html")


def viewResults(request, result_id):
    s1 = SequenceInput.objects.get(pk=result_id)
    sequence_name = s1.sequence_name
    sequence = s1.sequence
    temperature = s1.melt_Tm
    max_length = s1.max_Length
    oligo_concentration = s1.oligo_concentration
    monovalent_concentration = s1.monovalent_concentration
    output = OligoDesign.EquiTmOligo(sequence=sequence, oligo_conc=oligo_concentration, monovalent=monovalent_concentration, max_oligo_length=max_length, min_tm=temperature, min_oligo_length=5)
    output_oligos = output.main()
    if output_oligos[0] == True:
        output_oligos = output_oligos[1]
        output_oligo_profile = output.oligos_representation()
        primer_dimer = output.primer_dimer_check()
        secondary_structure = output.SecondaryStructureCheck()
        context = {"output_oligos": output_oligos, "result_id": result_id, "primer_dimer": primer_dimer, "secondary_structure": secondary_structure, "sequence_name": sequence_name, "result_id": result_id, "output_oligos_profile": output_oligo_profile}
        return render (request, "main/results.html", context)
    elif output_oligos[0] == False:
        context = {"output_oligos": output_oligos}
        return render(request, "main/error.html", context)
    else:
        return render(request, "404.html")
    try:
        sequence_output = SequenceOutput.objects.get(sequence__id=result_id)
    except ObjectDoesNotExist:
        s2 = SequenceOutput.objects.create(sequence=s1, output_oligos= output_oligo_profile, primer_dimers=primer_dimer, assembly_schem=output_oligos, secondary_structure=secondary_structure, all_results_url=reverse("main:results", args=[s1.id]))
    
        
def downloadResults(request, sequence_id):
    filename = "results.txt"
    results = SequenceOutput.objects.get(sequence__id=sequence_id)
    print(results)

    try:
        results = SequenceOutput.objects.get(sequence__id=sequence_id)
    except  ObjectDoesNotExist:
        output_oligos = ( False, "ERROR466: File Download Error", "An error occurred while creating your file. Please try again later.")
        context = {"output_oligos" : output_oligos}
        return render(request, "main/error.html", context) 
    except Exception as E:
        output_oligos = (False, "ERROR403: Forbidden Request", "You do not have access to the requested resource.")
        context = {"output_oligos" : output_oligos}
        return render(request, "main/error/html", context)
    
    oligos = results.output_oligos
    secondary_structure = results.secondary_structure
    content = assembly_scheme
    response = HttpResponse(content, content_type="text/plain")
    response['Content-Disposition'] = 'attachment; filename=results.txt'
    return response

