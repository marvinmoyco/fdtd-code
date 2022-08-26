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
import datetime


# Creating a class for the form submission...
class New_Simulation(forms.ModelForm):

    class Meta:
        model = Simulation
        fields = ['username',
                  'user_email',
                  'description',
                  'is_manualInput',
                  'input_csv',
                  'boundary_cond',
                  'excitation_method',
                  'custom_name',
                  'output_type',
                  'algorithm',
                  'multithreading',
                  'num_subdomains']


# Create your views here.
def index(request):
    return render(request, "web_interface/index.html")


def simulation(request,id):
    return render(request, "web_interface/simulation.html")


def add_simulation(request):
    
    #Check the submitted data in the POST request
    if request.method == "POST":
        datetime_str = datetime.datetime.now()
        print("Request.body")
        print(request.body)
        data = json.loads(request.body)
        print(data)
        print(f"username: {data['username']} | algo: {data['algorithm']}")
        # Create a new model based on the sent data

        #Store the data into the database
        return HttpResponseRedirect(reverse("web_interface:index"))
    # When the request is a GET request
    else:
        #return render(request, "web_interface/add_simulation.html")
        return JsonResponse({"error": "POST request required"},status=400)
    
    


    




