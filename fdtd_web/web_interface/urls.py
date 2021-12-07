from django.urls import path
from . import views

app_name='web_interface'
urlpatterns = [
    path("", views.index, name="index"),
    path("about",views.about, name="about"),
    
    #API Routes
    path("simulation/<int:id>",views.simulation,name="simulation"),
    path("simulation/add_simulation",views.add_simulation,name="add_simulation"),
]
