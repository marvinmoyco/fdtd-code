from django.db import models
from django.contrib.auth.models import AbstractUser
# Create your models here.



class User(AbstractUser):
    pass


class Simulation(models.Model):
    user = models.ForeignKey("User", on_delete=models.CASCADE, related_name="simulation")
    timestamp = models.DateTimeField(auto_now_add=True)
    