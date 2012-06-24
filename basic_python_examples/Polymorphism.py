#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-25
Description: Virtual Functions, polymorphism and Object Inheritance in Python

A simple example that shows how through inheritance from a base class and the
use of virtual functions we can create code that exhibits polymorphism.
That is, a common interface across different objects, such that the person
using our code may address different objects with a single interface.

In Python, all methods are virtual, so unlike in C++, there is no virtual
keyword needed to declare a method as being override-able.
'''

class Vehicle:
    def describe(self):
        print 'I am a Vehicle.'

class Automobile(Vehicle):
    def describe(self):
        print 'I am an Automobile.'

class Motorcycle(Vehicle):
    def describe(self):
        print 'I am a Motorcycle.'

class Truck(Vehicle):
    def describe(self):
        print 'I am a Truck.'
        
class GenericVehicle(Vehicle):
    pass

if __name__ == '__main__':
    vehicles = [Vehicle(), Automobile(), Motorcycle(), Truck(), GenericVehicle()]
    for vehicle in vehicles:
        vehicle.describe()

''' 
Python 2.7.3 (default, Apr 20 2012, 22:44:07) 
[GCC 4.6.3] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> execfile('Polymorphism.py')
I am a Vehicle.
I am an Automobile.
I am a Motorcycle.
I am a Truck.
I am a Vehicle.
'''