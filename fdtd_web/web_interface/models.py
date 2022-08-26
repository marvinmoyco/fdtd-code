from django.db import models
from django.contrib.auth.models import AbstractUser
# Create your models here.






class Simulation(models.Model):
    id = models.AutoField(primary_key=True)
    #Input data from the form
    username = models.CharField(max_length=255)
    user_email = models.EmailField(default=None)
    description = models.TextField(blank=True)
    #sim_param = models.JSONField()
    # sim_param contains the simulation parameters of the current simulation instance
    # sim_param contains the following parameters:
    # boundary condition: either dirichlet or pabc
    # source excitation method: either hard, soft, or tfsf method
    # custom output filename: filename for the output file (hdf5 and html file)
    # file type: csv, npy, or hdf5 (default)
    # algorithm: fdtd (default) or fdtd-schwarz
    # multithreading: true or false
    # number of subdomain: 0 to 64

    #Binary field to determine if the input parameter is by file upload or manual input
    is_manualInput = models.BooleanField(null=True,default=None)
    input_csv = models.FileField(null=True,default=None,upload_to='simulations/{id}')
    
    boundary_cond = models.CharField(null=True,default=None,max_length=255)
    excitation_method = models.CharField(null=True,default=None,max_length=255)
    custom_name = models.CharField(null=True,default=None,max_length=255)
    output_type = models.CharField(null=True,default=None,max_length=255)
    algorithm = models.CharField(null=True,default=None,max_length=255)
    multithreading = models.BooleanField(null=True,default=None)
    num_subdomains = models.PositiveIntegerField(default=1)
    sim_time = models.DurationField(default=0,null=True)
    #Simulation data (from the C++ file)
    
    timestamp = models.DateTimeField(auto_now_add=True)
    output_data = models.FileField(default=None,null=True,upload_to=f'simulations/{id}')
    log_data = models.FileField(null=True,default=None)

    def serialize(self):
        return {
             "id": self.id,
             "username": self.username,
             "user_email": self.user_email,
             "sim_description": self.description,
             "isManualInput": self.is_manualInput,
             "input_filepath": self.input_csv,
             #"sim_param": self.sim_param,
             "boundary_cond": self.boundary_cond,
             "excitation_method": self.excitation_method,
             "custom_name": self.custom_name,
             "output_type": self.output_type,
             "algorithm": self.algorithm,
             "multithreading": self.multithreading,
             "num_subdomains": self.num_subdomains,
             "output_filepath": self.output_data,
             "log_data": self.log_data,
             "timestamp": self.timestamp.strftime("%b %-d %Y, %-I:%M %p"),

        }
    