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
    if "sim_item" not in request.session:
        request.session["sim_item"] = None
    if "sim_list" not in request.session:
        request.session["sim_list"] = []
    return render(request, "web_interface/index.html")


def simulation(request):
    # Create a csv file based on the input file if the method is 'manual'
    input_filepath=""
    input_data = []
    if request.session["input_param"] != None:
        input_data = request.session["input_param"]
    else:
        return render(request, "web_interface/simulation.html", {
            "has_error" : True,
            "error_message" : "ERROR 001: No input data received from user! Please try again"
        })

    if input_data["input_type"] == 'manual':
        request.session["sim_item"]= loadJSONdata(request.session["input_param"])
        print(f"Sim object: {request.session['sim_item']}")
    else:
        # If the input model is already a file (csv format)
        pass

    # Check if a session variable of simulation object is created
    if request.session["sim_item"] == None:
        #return render(request, "web_interface/simulation.html", {
         #   "has_error" : True,
         #   "error_message" : "ERROR 002: No simulation object created during processing."
        #})
        return HttpResponseRedirect(reverse("web_interface:index"))
    #initialize computational domain
    request.session["sim_item"].init_comp_domain(spacer = 0,
                        inj_point = 0,
                        n_subdom = 1,
                        overlap = 5,
                        multithread = request.session["sim_item"].sim_param['multithreading_flag'],
                        algo=request.session["sim_item"].sim_param['algo'])
    return HttpResponseRedirect(reverse("web_interface:index"))


def add_simulation(request):
    
    #Check the submitted data in the POST request
    if request.method == "POST":
        datetime_str = datetime.datetime.now()
        print("Request.body")
        print(request.body)
        data = json.loads(request.body)
        request.session["input_param"] = data
        print(data)
        print(f"username: {data['username']} | algo: {data['algorithm']}")
        #Store the data into the database
        return HttpResponseRedirect(reverse("web_interface:simulation",))
     
    # When the request is a GET request
    else:
        # If the request is a GET request, load the HTML file for the form submission
        return render(request, "web_interface/add_simulation.html")
        
        #return JsonResponse({"error": "POST request required"},status=400)
    
    


    




