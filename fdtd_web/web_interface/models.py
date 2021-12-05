from django.db import models
from django.contrib.auth.models import AbstractUser
# Create your models here.






class Simulation(models.Model):
    id = models.AutoField(primary_key=True)
    user_email = models.EmailField(default=None)
    timestamp = models.DateTimeField(auto_now_add=True)
    output_data = models.FileField(default=None,upload_to=f'simulations/{id}')
    log_data = models.BinaryField(null=True,default=None)
    