import pygame
from pygame import gfxdraw


BG_COLOR = (40, 41, 35) # Dark Grey
color_R = (249, 36, 114) # RED
color_B = (103, 216, 239)
color_G = (166, 226, 43)
color_O = (253, 150, 34)
color_W = (248, 248, 239)
color_P = (172, 128, 255)
color_LG = (116, 112, 93)

class Circle:
    '''
    A function to give cursors and plots their shape and color.
    '''
    def __init__(self, x, y, color, radius):
        # Changed pos to contain both coordinates
        self.pos = (x, y)
        self.color = color
        self.radius = radius

    def radDist(self, plotCircle):
        return ((self.pos[0]-plotCircle.pos[0])**2 + (self.pos[1]-plotCircle.pos[1])**2)

    def draw_circle(self, surface):
        gfxdraw.aacircle(surface, self.pos[0], self.pos[1], self.radius, color_W)
        gfxdraw.filled_circle(surface, self.pos[0], self.pos[1], self.radius, color_W)
        gfxdraw.aacircle(surface, self.pos[0], self.pos[1], self.radius-2, self.color)
        gfxdraw.filled_circle(surface, self.pos[0], self.pos[1], self.radius-2, self.color)

    def draw_ring(self, surface, thickness):
        gfxdraw.aacircle(surface, self.pos[0], self.pos[1], self.radius+thickness//2, self.color)
        gfxdraw.filled_circle(surface, self.pos[0], self.pos[1], self.radius+thickness//2, self.color)
        gfxdraw.aacircle(surface, self.pos[0], self.pos[1], self.radius-thickness//2, BG_COLOR)
        gfxdraw.filled_circle(surface, self.pos[0], self.pos[1], self.radius-thickness//2, BG_COLOR)

class Button(pygame.Rect):
    '''
    Rectangular button graphics and functionality
    '''
    def __init__(self, x=0, y=0, width=100, height=20, text='Press', color=BG_COLOR):
        self.rectangle = pygame.Rect(x, y, width, height)

    def beingHovered(self, mouse_x, mouse_y):
        if (self.rectangle.left <= mouse_x <= self.rectangle.right 
            and self.rectangle.top <= mouse_y <= self.rectangle.bottom):
            return True
        else:
            return False

