#!/usr/bin/env python

import pygame
import random
import sys
import Vector


class BoardObject(object):
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        self.xmax = xmax
        self.ymax = ymax
        self.xpos = xpos
        self.ypos = ypos
        self.heading = Vector.Vector(speed, bearing)
        self.orientation = 0.0
        self.size = 10.0

    def update(self):
        '''
        Given the current (x,y) position and the heading, calculate the next
        (x,y) position.
        '''
        self.xpos = (self.xpos + self.heading.x) % self.xmax
        self.ypos = (self.ypos + self.heading.y) % self.ymax

    def draw(self):
        return pygame.draw.polygon

class Ship(BoardObject):
    THRUST = 0.5
    ROTATE = 10
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed, bearing)

    def thrust(self):
        self.heading = self.heading + Vector.Vector(Ship.THRUST, self.orientation)

    def reversethrust(self):
        self.heading = self.heading + Vector.Vector(-Ship.THRUST, self.orientation)

    def rotate_left(self):
        self.orientation = (self.orientation - Ship.ROTATE) % 360.0

    def rotate_right(self):
        self.orientation = (self.orientation + Ship.ROTATE) % 360.0

    def points(self):
        '''Based on the orientation, return the point list to properly orient the ship in space.'''
        front      = Vector.Vector(10, self.orientation)
        left_rear  = Vector.Vector(10, (self.orientation-120) % 360)
        right_rear = Vector.Vector(10, (self.orientation+120) % 360)
        rear       = Vector.Vector(-1, (self.orientation+180) % 360)
        return (((self.xpos+front.x), (self.ypos+front.y)),
                ((self.xpos+left_rear.x), (self.ypos+left_rear.y)),
                ((self.xpos+rear.x), (self.ypos+rear.y)),
                ((self.xpos+right_rear.x), (self.ypos+right_rear.y)),)

if __name__ == '__main__':
    pygame.init()
    fpsClock = pygame.time.Clock()
    pygame.display.set_caption('Asteroids')
    BLACK = pygame.Color(0,0,0)
    WHITE = pygame.Color(255,255,255)
    xmax = 640
    ymax = 480
    screen = pygame.display.set_mode((xmax,ymax))
    xcenter = xmax/2
    ycenter = ymax/2

    myShip = Ship(xmax, ymax, xcenter, ycenter, 0.0, 0.0)
    Board = list()
    Board.append(myShip)

    while True:
        screen.fill(BLACK)

        for i in Board:
            i.update()
            i.draw()(screen, WHITE, i.points(), 1)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_LEFT:
                    myShip.rotate_left()
                elif event.key == pygame.K_RIGHT:
                    myShip.rotate_right()
                elif event.key == pygame.K_UP:
                    myShip.thrust()
                elif event.key == pygame.K_DOWN:
                    myShip.reversethrust()
        pygame.display.update()
        fpsClock.tick(60)
