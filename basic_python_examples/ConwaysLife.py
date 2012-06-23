#!/usr/bin/env python

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
'''

import random
import sys

class Life:
    def __init__(self, nrow=4, ncol=4):
        self.nrow = nrow
        self.ncol = ncol
        self.size = self.nrow * self.ncol
        self.current_grid = [random.randint(0,1) for i in xrange(self.size)]
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
        self.current_grid = self.next_grid
        self.generation += 1
    
    def show_state(self):
        index = 0
        char = '_'
        for cell in self.current_grid:
            if self.current_grid[index] == 1:
                char = '+'
            else:
                char = '_'
            if index % self.ncol == (self.ncol - 1):
                sys.stdout.write('{0}\n'.format(char))
            else:
                sys.stdout.write('{0}'.format(char))
            index += 1