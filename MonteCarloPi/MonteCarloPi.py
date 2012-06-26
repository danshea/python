#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-26
Description: Graphical Monte Carlo Simulation to calculate pi

If we inscribe a square with a circle and place points randomly within the square,
we can see that four times the ratio of points that fall within the circle in relation to the
total number of points will approximate pi.

Area_square = 4r**2
Area_circle = pi*r**2

pi*r**2 / 4r**2 = num_points_in_circle / total_points

pi = 4*num_points_circle / total_points
'''

import math
import random
import wx

class MonteCarloPi(wx.Frame):
    ID_TIMER = 1
    speed = 10
    
    def __init__(self, parent, title, size=500):   
        super(MonteCarloPi, self).__init__(parent, title = title, size = (size, size))
        self.timer = wx.Timer(self, MonteCarloPi.ID_TIMER)
        
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_TIMER, self.OnTimer, id=MonteCarloPi.ID_TIMER)

        self.size = size
        self.radius = size/2
        self.points_in_circle = list()
        self.points_outside_circle = list()
        self.pi = 'pi={0}'.format(0.0)
        
        self.start()

        self.Centre()
        self.Show()
    
    def OnTimer(self, event):
        if event.GetId() == MonteCarloPi.ID_TIMER:
            self.next()
            self.Refresh()
        else:
            event.Skip()
    
    def start(self):
        self.isStarted = True
        self.timer.Start(MonteCarloPi.speed)
    
    def next(self):
        self.generate_point()
        self.estimate_pi()
    
    def in_circle(self, x, y):
        if math.sqrt((self.radius - x)**2 + (self.radius - y)**2) < self.radius:
            return True
        else:
            return False
        
    def estimate_pi(self):
         self.pi = 'pi = {0}'.format(4.0 * len(self.points_in_circle) / (len(self.points_in_circle) + len(self.points_outside_circle)))
        
    def generate_point(self):
        x, y = random.randint(0,self.size), random.randint(0,self.size) 
        if self.in_circle(x, y):
            self.points_in_circle.append((x,y))
        else:
            self.points_outside_circle.append((x, y))
    
    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        
        dc.SetPen(wx.Pen('#000000'))
        dc.SetBrush(wx.Brush('#0C0CF0'))
        dc.DrawRectangle(0,0,self.size,self.size)
        
        dc.SetBrush(wx.Brush('#7EE7F2'))
        dc.DrawCircle(self.radius, self.radius, self.radius)
        
        for point in self.points_in_circle:
            dc.SetPen(wx.Pen('#000000'))
            dc.DrawPoint(*point)
        
        for point in self.points_outside_circle:
            dc.SetPen(wx.Pen('#FCCACA'))
            dc.DrawPoint(*point)
            
        dc.SetTextForeground('#FFFFFF')
        dc.DrawText(self.pi, 10, 10)

if __name__ == "__main__":
    app = wx.App()
    MonteCarloPi = MonteCarloPi(None, title = 'Monte Carlo Simulation of Pi')
    app.MainLoop()