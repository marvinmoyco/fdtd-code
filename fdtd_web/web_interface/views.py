from django.contrib.auth import authenticate, login, logout
from django.db import IntegrityError
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django.urls import reverse
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.decorators import login_required
import json
from django import forms
import os
from .models import *

# Creating a class for the form submission...
class New_Simulation(forms.Form):
    username = forms.CharField()
    description = forms.CharField()
    
    #Selects whether the method is by uploading a csv file or manual entry
    selected_method = forms.ChoiceField()
    #csv_filepath = forms.FilePathField()

    fmax = forms.FloatField()
    source_type = forms.ChoiceField()
    num_layers = forms.IntegerField()

    
    #Get the values based on the num_layers
    #layer_size = [forms.FloatField() for x in range(num_layers)]
    ##mu = [forms.FloatField() for x in range(num_layers)]
    #epsilon = [forms.FloatField() for x in range(num_layers)]

    boundary_cond = forms.ChoiceField()
    excitation_method = forms.ChoiceField()
    custom_name = forms.CharField()
    output_file = forms.ChoiceField()
    algorithm = forms.ChoiceField()
    multithreading = forms.BooleanField()
    num_subdomains = forms.IntegerField()

# Create your views here.
def index(request):
    return render(request, "web_interface/index.html")


def simulation(request,id):
    return render(request, "web_interface:simulation.html")

@csrf_exempt
def add_simulation(request):
    
    #Check the submitted data in the POST request
    if request.method != "POST":
        return JsonResponse({"error": "POST request required"},status=400)

        
    data = json.loads(request.body)
    print(data)
    if data != None:
        print("HELOOOOOOOOOOOOOOOOOOOO")
    else:
        return render(request, "web_interface/index.html", {
            "has_error": True,
            "error_message": "Submitted data is either incomplete/invalid"
        })
        

    




