#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: djs
Date: 2012-06-24
Description: Conway's Game of Life

The universe of the Game of Life is an infinite two-dimensional orthogonal grid of square cells,
each of which is in one of two possible states, alive or dead. Every cell interacts with its
eight neighbours, which are the cells that are horizontally, vertically, or diagonally adjacent.
At each step in time, the following transitions occur:

    Any live cell with fewer than two live neighbours dies, as if caused by under-population.
    Any live cell with two or three live neighbours lives on to the next generation.
    Any live cell with more than three live neighbours dies, as if by overcrowding.
    Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.

The initial pattern constitutes the seed of the system. The first generation is created by
applying the above rules simultaneously to every cell in the seed-births and deaths occur
simultaneously, and the discrete moment at which this happens is sometimes called a tick
(in other words, each generation is a pure function of the preceding one). The rules continue
to be applied repeatedly to create further generations.

This version of the implmentation cheats the infinite two-dimensional orthogonal grid requirement
by implementing a toroidal array.  The left/right edges wrap and the top/bottom edges wrap so that
cells on the boundaries of the array reference cells on the opposing edges as neighbors.

The gui is implmented in wxPython for display of the simulation.
'''

import os
import os.path
import time
import sys
import random
import wx

class Life():
    def __init__(self, nrow=50, ncol=50, initial_grid = None):
        self.nrow = nrow
        self.ncol = ncol
        self.size = self.nrow * self.ncol
        if initial_grid == None:
            self.current_grid = [random.randint(0,1) for i in xrange(self.size)]
        else:
            self.current_grid = initial_grid
        self.next_grid = [0 for i in xrange(self.size)]
        self.generation = 0
    
    def _up(self, index):
        return (index - self.ncol) % self.size
    
    def _down(self, index):
        return (index + self.ncol) % self.size
    
    def _left(self, index):
        if index % self.ncol == 0:
            return index + (self.ncol - 1)
        else:
            return index - 1
    
    def _right(self, index):
        if index % self.ncol == (self.ncol - 1):
            return index - (self.ncol - 1)
        else:
            return index + 1
    
    def _lup(self, index):
        return self._left(self._up(index))
    
    def _rup(self, index):
        return self._right(self._up(index))
        
    def _ldown(self, index):
        return self._left(self._down(index))
    
    def _rdown(self, index):
        return self._right(self._down(index))
    
    def _get_neighbors(self, index):
        return [self.current_grid[self._lup(index)], self.current_grid[self._up(index)], self.current_grid[self._rup(index)],
                self.current_grid[self._left(index)], self.current_grid[self._right(index)],
                self.current_grid[self._ldown(index)], self.current_grid[self._down(index)], self.current_grid[self._rdown(index)]]
    
    def next(self):
        index = 0
        for i in self.current_grid:
            neighbors = self._get_neighbors(index)
            if self.current_grid[index] == 1:
                # Any live cell with fewer than two live neighbours dies, as if caused by under-population.
                if sum(neighbors) < 2:
                    self.next_grid[index] = 0
                # Any live cell with two or three live neighbours lives on to the next generation.
                elif 2 <= sum(neighbors) <= 3:
                    self.next_grid[index] = 1
                # Any live cell with more than three live neighbours dies, as if by overcrowding.
                elif sum(neighbors) > 3:
                    self.next_grid[index] = 0
            else:
                # Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.
                if sum(neighbors) == 3:
                    self.next_grid[index] = 1
            index += 1
        # mistakenly was pointing the current_grid to the next grid instead of making a copy
        self.current_grid = self.next_grid[:]
        self.generation += 1
    
    def show_state(self):
        index = 0
        char = ''
        for cell in self.current_grid:
            if self.current_grid[index] == 1:
                char = '+'
            else:
                char = '-'
            if index % self.ncol == (self.ncol - 1):
                sys.stdout.write('{0}\n'.format(char))
            else:
                sys.stdout.write('{0}'.format(char))
            index += 1

class LifeApp(wx.Frame):
    ID_TIMER = 1
    speed = 1000
    
    def __init__(self, parent, title, xsize=500, ysize=500, cell_size=10, initial_grid = None):   
        super(LifeApp, self).__init__(parent, title = title, size = (xsize, ysize))
        self.timer = wx.Timer(self, LifeApp.ID_TIMER)
        
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_TIMER, self.OnTimer, id=LifeApp.ID_TIMER)

        self.cell_size = cell_size
        self.LifeSim = Life(xsize/self.cell_size, ysize/self.cell_size, initial_grid)

        self.start()

        self.Centre()
        self.Show()
    
    def OnTimer(self, event):
        if event.GetId() == LifeApp.ID_TIMER:
            self.LifeSim.next()
            self.Refresh()
        else:
            event.Skip()
    
    def start(self):
        self.isStarted = True
        self.timer.Start(LifeApp.speed)
    
    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        index = 0
        for cell in self.LifeSim.current_grid:
            i = (index / self.LifeSim.ncol)
            j = (index % self.LifeSim.ncol)
            dc.SetPen(wx.Pen((0, 0, 0)))
            if self.LifeSim.current_grid[index] == 0:
                dc.SetBrush(wx.Brush((0, 0, 0)))
            else:
                dc.SetBrush(wx.Brush((35, 142, 35)))
            dc.DrawRectangle(self.cell_size*j, self.cell_size*i, self.cell_size, self.cell_size)
            index += 1

def usage():
    print 'usage: {0} [filename]'.format(sys.argv[0])
    sys.exit(0)

if __name__ == "__main__":
    app = wx.App()
    if len(sys.argv) == 2:
        if os.path.isfile(sys.argv[1]):
            with open(sys.argv[1], 'r') as fh:
                grid = ''
                for line in fh:
                    grid += line.strip()
                grid = eval(grid)
            if type(grid) is list:
                lifeapp = LifeApp(None, title = 'Conway\'s Life', initial_grid = grid)
            else:
                usage()
    elif len(sys.argv) == 1:
        lifeapp = LifeApp(None, title = 'Conway\'s Life')
    else:
        usage()
    app.MainLoop()    
