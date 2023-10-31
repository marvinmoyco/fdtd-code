from django.urls import path
from django.conf.urls.static import static
import sys
from . import views

app_name='web_interface'
urlpatterns = [
    path("", views.index, name="index"),
    
    #API Routes
    path("simulation",views.simulation,name="simulation"),
    path("new",views.add_simulation,name="add_simulation"),
]
