#!/usr/bin/env python3
#basic_classes.py

class Circle():
    def __init__(self, color, radius):
        self.color = color
        self.radius = radius
        
    def diameter(self):
        '''finds the diameter of a circle'''
        return 2 * self.radius
    
    def circumference(self):
        ''' finds the circumference of a circle'''
        return 3.14 * 2 * self.radius
        
    def isRed(self):
        if self.color == "red":
            print("True")
        else:
            return False
wheel = Circle("red", 14)
wheel.diameter()

class GraduateStudent():
    def __init__(self, first_name, last_name, year, major):
        self.first_name = first_name
        self.last_name = last_name
        self.year = year
        self.major = major

    def year_matriculated(self):
        ''' finds year student was matriculted'''
        return 2020 - self.year

person_1 = GraduateStudent("Rebecca", "Meyer", 2, "computer science")
person_1.year_matriculated()

person_2 = GraduateStudent("Stefan", "Kaluziak","5","Marine and environmental science")
person_2.year_matriculated()
