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
        # 4 points of the hitbox for collision detection
        # top left, bottom left, top right, bottom right
        self.hitbox = (((self.xpos-self.size)%self.xmax, (self.ypos-self.size)%self.ymax),
                       ((self.xpos-self.size)%self.xmax, (self.ypos+self.size)%self.ymax),
                       ((self.xpos+self.size)%self.xmax, (self.ypos-self.size)%self.ymax),
                       ((self.xpos+self.size)%self.xmax, (self.ypos+self.size)%self.ymax),)

    def update(self):
        '''
        Given the current (x,y) position and the heading, calculate the next
        (x,y) position.  Update the hitbox coordinates.
        '''
        self.xpos = (self.xpos + self.heading.x) % self.xmax
        self.ypos = (self.ypos + self.heading.y) % self.ymax
        self.hitbox = (((self.xpos-self.size)%self.xmax, (self.ypos-self.size)%self.ymax),
                       ((self.xpos-self.size)%self.xmax, (self.ypos+self.size)%self.ymax),
                       ((self.xpos+self.size)%self.xmax, (self.ypos-self.size)%self.ymax),
                       ((self.xpos+self.size)%self.xmax, (self.ypos+self.size)%self.ymax),)

    def draw(self):
        return pygame.draw.polygon

class Bullet(BoardObject):
    SPEED = 5.0
    TTL = 75
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed+Bullet.SPEED, bearing)
        self.ttl = Bullet.TTL
        self.orientation = bearing
        self.size = 2

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
    SPAWN_SIZE = 40
    SCORE = 100
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed, bearing)
        self.size = Asteroid.SPAWN_SIZE
        self.orientation = random.randint(0,359)
        self.score = Asteroid.SCORE

    def points(self):
        p1    = Vector.Vector(self.size, self.orientation)
        p2    = Vector.Vector(self.size, (self.orientation+72.5) % 360)
        p3    = Vector.Vector(self.size, (self.orientation+145) % 360)
        p4    = Vector.Vector(self.size, (self.orientation+217.5) % 360)
        p5    = Vector.Vector(self.size, (self.orientation+290) % 360)
        return (((self.xpos+p1.x), (self.ypos+p1.y)),
                ((self.xpos+p2.x), (self.ypos+p2.y)),
                ((self.xpos+p3.x), (self.ypos+p3.y)),
                ((self.xpos+p4.x), (self.ypos+p4.y)),
                ((self.xpos+p5.x), (self.ypos+p5.y)),)

class Ship(BoardObject):
    THRUST = 0.5
    ROTATE = 5
    def __init__(self, xmax, ymax, xpos, ypos, speed, bearing):
        BoardObject.__init__(self, xmax, ymax, xpos, ypos, speed, bearing)
        self.weapon = Bullet
        self.lives = 2

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

    def respawn(self):
        self.xpos = self.xmax / 2
        self.ypos = self.ymax / 2
        self.orientation = 0.0
        self.heading = Vector.Vector(0, 0)

class Board(object):

    def __init__(self):
        self.board = list()
        self.level = 1
        self.score = 0

    def game_over(self):
        '''game is over'''
        pass

    def collision_check(self, A, B):
        '''collision_check - if the hitboxes of A and B overlap, return True as a collision has occurred.'''
        for point in A.hitbox:
            if B.hitbox[0][0] <= point[0] <= B.hitbox[3][0]:
                if B.hitbox[0][1] <= point[1] <= B.hitbox[3][1]:
                    return True
        return False

    def detect_collisions(self):
        to_remove = set()
        to_add = list()
        for i in self.board:
            for j in self.board:
                if self.collision_check(i, j) == True:
                    if (type(i) is Ship and type(j) is Asteroid) or (type(i) is Asteroid and type(j) is Ship):
                        if type(i) is Ship:
                            i.lives -= 1
                            if i.lives < 0:
                                self.game_over()
                            else:
                                i.respawn()
                        else:
                            j.lives -= 1
                            if j.lives < 0:
                                self.game_over()
                            else:
                                j.respawn()
                    elif (type(i) is Bullet and type(j) is Asteroid) or (type(i) is Bullet and type(j) is Asteroid):
                        to_remove.add(i)
                        to_remove.add(j)
                        new_heading = i.heading + j.heading
                        if type(i) is Asteroid:
                            self.score += i.score
                            new_heading.r = i.heading.r
                            if i.size / 2 > 2:
                                lchild = Asteroid(i.xmax, i.ymax, i.xpos, i.ypos, new_heading.r, new_heading.theta)
                                rchild = Asteroid(i.xmax, i.ymax, i.xpos, i.ypos, new_heading.r, -new_heading.theta)
                                lchild.score = i.score * 2
                                rchild.score = i.score * 2
                                lchild.size = i.size / 2
                                rchild.size = i.size / 2
                                speed_split = (random.randint(5,95) / 100.0)
                                lchild.heading.r = i.heading.r * speed_split
                                rchild.heading.r = i.heading.r * (1 - speed_split)
                                to_add.append(lchild)
                                to_add.append(rchild)
                        elif j.size/2 > 2:
                            self.score += j.score
                            lchild = Asteroid(j.xmax, j.ymax, j.xpos, j.ypos, new_heading.r, new_heading.theta)
                            rchild = Asteroid(j.xmax, j.ymax, j.xpos, j.ypos, new_heading.r, -new_heading.theta)
                            lchild.score = j.score * 2
                            rchild.score = j.score * 2
                            lchild.size = j.size / 2
                            rchild.size = j.size / 2
                            speed_split = (random.randint(5,95) / 100.0)
                            lchild.heading.r = j.heading.r * speed_split
                            rchild.heading.r = j.heading.r * (1 - speed_split)
                            to_add.append(lchild)
                            to_add.append(rchild)
        for i in to_remove:
            self.board.remove(i)
        for i in to_add:
            self.board.append(i)

if __name__ == '__main__':
    pygame.init()
    fpsClock = pygame.time.Clock()
    pygame.display.set_caption('Asteroids')
    pygame.key.set_repeat(150,50)
    BLACK = pygame.Color(0,0,0)
    WHITE = pygame.Color(255,255,255)
    RED = pygame.Color(255,0,0)
    xmax = 800
    ymax = 600
    screen = pygame.display.set_mode((xmax,ymax))
    xcenter = xmax/2
    ycenter = ymax/2

    myShip = Ship(xmax, ymax, xcenter, ycenter, 0.0, 0.0)
    board = Board()
    board.board.append(myShip)
    for i in xrange(board.level):
        board.board.append(Asteroid(xmax, ymax, random.randint(0, xmax/2 - xmax/4), random.randint(0, ymax/2 - ymax/4), random.randint(1,5)/2.0, random.randint(0,359)))
        board.board.append(Asteroid(xmax, ymax, random.randint(0, xmax/2 - xmax/4), random.randint(ymax/2 + ymax/4, ymax), random.randint(1,5)/2.0, random.randint(0,359)))
        board.board.append(Asteroid(xmax, ymax, random.randint(xmax/2 + xmax/4, xmax), random.randint(0, ymax/2 - ymax/4), random.randint(1,5)/2.0, random.randint(0,359)))
        board.board.append(Asteroid(xmax, ymax, random.randint(xmax/2 + xmax/4, xmax), random.randint(ymax/2 + ymax/4, ymax), random.randint(1,5)/2.0, random.randint(0,359)))

    while True:
        screen.fill(BLACK)

        for i in board.board:
            # update positions on the board
            i.update()

        # perform collision detection
        board.detect_collisions()

        # redraw the board
        index = 0
        for i in board.board:
            if type(i) is Bullet:
                if i.ttl > 0:
                    i.draw()(screen, RED, i.points(), 1)
                else:
                    board.board.pop(index)
            else:
                i.draw()(screen, WHITE, i.points(), 1)
            index += 1

        # capture keyboard events
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
                    board.board.append(myShip.shoot())
        pygame.display.update()
        fpsClock.tick(60)
