from django.urls import path
from . import views

app_name='web_interface'
urlpatterns = [
    path("", views.index, name="index"),
    
    #API Routes
    path("simulation/<int:id>",views.simulation,name="simulation"),
    path("new",views.add_simulation,name="add_simulation"),
]
