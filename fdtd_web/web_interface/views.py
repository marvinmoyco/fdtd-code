from django.contrib.auth import authenticate, login, logout
from django.db import IntegrityError
from django.http import HttpResponse, HttpResponseRedirect, Http404, HttpResponseForbidden
from django.shortcuts import render
from django.urls import reverse
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.decorators import login_required
import json
from django import forms
import os
from .models import *
import csv
import sys
from .simulation import *
from .utility import loadJSONdata
import datetime


# fix windows registry stuff
import mimetypes
mimetypes.add_type('application/javascript', '.js')
mimetypes.add_type('text/css', '.css')


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

        # Create a csv file based on the input file if the method is 'manual'
        input_filepath=""
        if data['input_type'] == 'manual':
            sim = loadJSONdata(data)
            print(f'layer_size: {data["layer_size"]}')
            print(f"mu: {data['mu']}")
            print(f"epsilon: {data['epsilon']}")
            
            




        else:
            # If the input model is already a file (csv format)
            pass

        # Create a new model based on the sent data
        


        #Store the data into the database
        #return HttpResponseRedirect(reverse("web_interface:simulation", args=(simulation.id,)))
        return render(request, "web_interface/simulation.html", {"sim_id": 1})
    # When the request is a GET request
    else:
        #return render(request, "web_interface/add_simulation.html")
        return HttpResponseForbidden()
        #return JsonResponse({"error": "POST request required"},status=400)
    
    


    




