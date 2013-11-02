#!/usr/bin/env python3

class Node(object):
    def __init__(self, xmax, ymax, index):
        self.x = index % xmax
        self.y = index // ymax
        self.index = index
        self.walkable = True
        self.cost = 1.0

    def up(self):
        y = (self.y - 1) % self.ymax
        x = self.x
        return (self.xmax * y) + x

    def down(self):
        y = (self.y + 1) % self.ymax
        x = self.x
        return (self.xmax * y) + x

    def left(self):
        y = self.y
        x = (self.x - 1) % self.xmax
        return (self.xmax * y) + x

    def right(self):
        y = self.y
        x = (self.x + 1) % self.xmax
        return (self.xmax * y) + x

    def up_left(self):
        y = (self.y - 1) % self.ymax
        x = (self.x - 1) % self.xmax
        return (self.xmax * y) + x

    def down_left(self):
        y = (self.y + 1) % self.ymax
        x = (self.x - 1) % self.xmax
        return (self.xmax * y) + x

    def up_right(self):
        y = (self.y - 1) % self.ymax
        x = (self.x + 1) % self.xmax
        return (self.xmax * y) + x

    def down_rightt(self):
        y = (self.y + 1) % self.ymax
        x = (self.x + 1) % self.xmax
        return (self.xmax * y) + x

    def neighbors(self):
        return [self.up(), self.up_right(), self.right(),
                self.down_right(), self.down(), self.down_left(),
                self.left(), self.up_left(),]
    
