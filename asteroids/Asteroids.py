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

class Bullet(BoardObject):
    SPEED = 5.0
    TTL = 75
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed+Bullet.SPEED, bearing)
        self.ttl = Bullet.TTL
        self.orientation = bearing

    def update(self):
        BoardObject.update(self)
        self.ttl -= 1

    def points(self):
        front = Vector.Vector(1, self.orientation)
        lr    = Vector.Vector(1, (self.orientation-120) % 360)
        rr    = Vector.Vector(1, (self.orientation+120) % 360)
        return (((self.xpos+front.x), (self.ypos+front.y)),
                ((self.xpos+lr.x), (self.ypos+lr.y)),
                ((self.xpos+rr.x), (self.ypos+rr.y)),)

class Asteroid(BoardObject):
    SPEED = 1
    SPAWN_SIZE = 50
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed, bearing)
    

class Ship(BoardObject):
    THRUST = 0.5
    ROTATE = 10
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed, bearing)
        self.weapon = Bullet

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

    def shoot(self):
        speed = Vector.Vector(0, self.orientation) + self.heading
        return self.weapon(self.xmax, self.ymax, self.xpos, self.ypos, speed.r, self.orientation)

if __name__ == '__main__':
    pygame.init()
    fpsClock = pygame.time.Clock()
    pygame.display.set_caption('Asteroids')
    BLACK = pygame.Color(0,0,0)
    WHITE = pygame.Color(255,255,255)
    RED = pygame.Color(255,0,0)
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

        index = 0
        for i in Board:
            i.update()
            if type(i) is Bullet:
                if i.ttl > 0:
                    i.draw()(screen, RED, i.points(), 1)
                else:
                    Board.pop(index)
            else:
                i.draw()(screen, WHITE, i.points(), 1)
            index += 1

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
                elif event.key == pygame.K_SPACE:
                    Board.append(myShip.shoot())
        pygame.display.update()
        fpsClock.tick(60)
