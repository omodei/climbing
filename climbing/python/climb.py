#!/usr/bin/env python

import tkinter as tkr
import time
from matplotlib import pyplot as plt
import numpy as np
from data import RopeData
#Rope model from https://www.sigmadewe.com/fileadmin/user_upload/pdf-Dateien/SEILPHYSIK.pdf

class Geometry():
    def __init__(self,W=20,H=50,scale=1):
        self.H = H
        self.W = W
        self.m2pix = scale * 20
        # world coordinates:
        self.center=(W/2.,H/2.)
        self.top    =  H / 2.
        self.bottom = -H / 2.
        self.left   = -W / 2.
        self.right  =  W / 2.

        # create canvas:
        self.tk = tkr.Tk()
        self.canvas = tkr.Canvas(self.tk,width=W * self.m2pix,height=H * self.m2pix,bg='white')
        self.canvas.grid()
        self.text = None

    def DisplayTime(self,text):
        if self.text: self.canvas.delete(self.text)
        self.text = self.canvas.create_text(0.9 * self.W * self.m2pix, 0.1 * self.H * self.m2pix, text=text)


    def world2coord(self, pos):
        '''
        Transorms world coordinates in meters to pixels
        x=0,y=0 is the center of the canvas. +y are toward the top, -y are bottom.
        :param pos: position tuple (x,y) in meters
        :return: tupla of position x,y in pixels
        '''
        return np.array([pos[0] + self.center[0], self.center[1] - pos[1]]) * self.m2pix

    def coord2world(self, pos):
        '''
        Transorms pixels coordinates in worls coordinate in meters
        x=0,y=0 is the center of the canvas. +y are toward the top, -y are bottom.
        :param pos: position tuple (x,y) in pixels
        :return: tupla of position x,y in meters
        '''
        return np.array([pos[0] / self.m2pix - self.center[0], self.center[1] - pos[1] / self.m2pix])

    def DrawLine(self,x0,y0,x1,y1,**opt):
        pos = list(np.append(self.world2coord((x0,y0)), self.world2coord((x1,y1))))
        self.canvas.create_line(pos, opt)

    def DrawReferenceGrid(self,step=10):
        '''
        Drwa
        :param step:
        :return:
        '''
        ys=np.arange(self.bottom,self.top,step)
        for y in ys:
            self.DrawLine(self.left,  y, self.right, y, fill='gray')
            pass
        self.DrawLine(self.left, 0, self.right, 0, fill='red')

    def DrawQuickDraw(self,x0,y0,angle=0, **opt):
        '''
        Draw a quickdraw
        :param x0: x position of center of the carabiner with rope
        :param y0: y position of center of the carabiner with rope
        :param angle: angle with respect the vertical.
        :param opt: optional argument, for example fill='red'
        :return:
        '''
        #image = Image.open('images/QD.png')
        dx     = 0.04
        bone   = .15
        angle *=np.pi/180.

        #first carabiner(rope):
        pos1  = self.world2coord((x0 - dx, y0 - dx))
        pos2  = self.world2coord((x0 + dx, y0 + dx))
        pos  = list(np.append(pos1,pos2))
        self.canvas.create_oval(pos, opt)
        #bone

        pos1  = self.world2coord((x0, y0))
        pos2  = self.world2coord((x0 + bone * np.sin(angle), y0 + bone * np.cos(angle)))
        pos   = list(np.append(pos1,pos2))
        self.canvas.create_line(pos, opt)

        #seond carabiner (bolt):
        pos1  = self.world2coord((x0 + bone * np.sin(angle) - dx, y0 + bone * np.cos(angle) - dx ))
        pos2  = self.world2coord((x0 + bone * np.sin(angle) + dx, y0 + bone * np.cos(angle) + dx))
        pos   = list(np.append(pos1,pos2))
        self.canvas.create_oval(pos, opt)




    def mainloop(self):
        self.tk.mainloop()

class Climber():
    def __init__(self,geometry,mass,x,y):
        self.geo = geometry
        self.m   = mass
        self.dx  = 0.25
        pos1 = self.geo.world2coord((x - self.dx, y - self.dx))
        pos2 = self.geo.world2coord((x + self.dx, y + self.dx))
        pos = list(np.append(pos1, pos2))

        self.climber = self.geo.canvas.create_oval(pos, fill='blue')
        self.pos   = self.get_position()
        self.speed = np.array([0. , 0.])
        self.accel = np.array([0. , 0.]) # gravity add by default
        self.force_vector = None
        pass

    def set_force(self,F):
        self.accel=F/self.m

    def set_speed(self,vx,vy):
        self.speed = np.array([vx, vy])

    def set_accel(self,ax,ay):
        self.accel = np.array([ax, ay])

    def get_force(self):
        return self.m * self.accel

    def get_speed(self):
        return self.speed

    def get_force_magnitude(self):
        return np.linalg.norm(self.get_force())

    def get_force_direction(self):
        return self.get_force()/np.linalg.norm(self.get_force())


    def get_canvas_position(self):
        pos =  self.geo.canvas.coords(self.climber)

        x   = 0.5*(pos[0] + pos[2])
        y   = 0.5*(pos[1] + pos[3])
        return x,y

    def get_position(self):
        x,y = self.get_canvas_position()
        return self.geo.coord2world((x,y))

    def update(self,dt):
        v=self.speed
        a=self.accel
        v+= a * dt
        self.set_speed(v[0],v[1])
        self.geo.canvas.move(self.climber, v[0] * self.geo.m2pix * dt, -v[1] * self.geo.m2pix * dt)
        self.pos   = self.get_position()

    def draw_force(self,scale,refresh=True):
        x,y   = self.get_canvas_position()
        dx,dy = self.get_force() / scale
        if refresh and self.force_vector: self.geo.canvas.delete(self.force_vector)
        self.force_vector = self.geo.canvas.create_line(x,y,
                                                        x + dx,
                                                        y - dy, fill='red',arrow=tkr.LAST,width=3)


class Rope():
    def __init__(self,geometry):
        self.geo=geometry
        self.rope_color='green'
        self.rope      = None
        self.rope_points = []
        #self.k1 = 5916.8
        #self.k2 = 1620.0
        #self.c =   1376.0

        self.k1 = 6000
        self.k2 = 1500
        self.c  = 1500

        self.force_vector = None
        self.tension      = None
        self.tension_direction = None

    def set_constants(self,k1,k2,c):
        self.k1 = k1
        self.k2 = k2
        self.c  = c

    def add_anchor(self,x,y,angle):
        self.anchor_pos_world = np.array([x,y])
        self.anchor_pos       = self.geo.world2coord(self.anchor_pos_world)
        self.rope_points.append(self.anchor_pos_world )
        self.geo.DrawQuickDraw(x, y, angle=angle, fill='red')


    def add_quickdraw(self,x,y,angle):
        self.rope_points.append(np.array([x,y]))
        self.geo.DrawQuickDraw(x,y, angle=angle, fill='black')


    def add_climber(self,climber):
        self.climber = climber
        self.rope_points.append(self.climber.get_position())

    def compute_length(self):
        d=0
        for a,b in zip(self.rope_points[:-1],self.rope_points[1:]):
            d += np.linalg.norm(np.array(a)-np.array(b))
        return d

    def compute_tension_direction(self):
        '''
        Compute the projection on the x and y axis
        :return: cos(alpha),sin(alpha) with alpha the angle with the vertical
        '''
        a = np.array(self.rope_points[-1])
        b = np.array(self.rope_points[-2])
        d=np.linalg.norm(a-b)
        return (b-a)/d

    #def compute_fall_factor(self):
    #    return self.climber.get_position()/


    def set_rope_length(self,length):
        self.length = length


    def update(self):
        self.rope_points[-1] = self.climber.get_position()
        self.get_stretch()
        self.compute_tension()

    def draw(self,scale=1.0,refresh=True):

        if self.rope: self.geo.canvas.delete(self.rope)
        if not self.stretch: color = self.rope_color
        else: color = 'red'
        self.rope_points_canvas = []
        for xy in self.rope_points:
            self.rope_points_canvas.append(self.geo.world2coord((xy[0],xy[1])))

        list_of_points = list(np.array(self.rope_points_canvas).reshape(len(self.rope_points_canvas)*2))

        self.rope = self.geo.canvas.create_line(list_of_points,fill=color,width=2)

        if scale:
            if refresh and self.force_vector: self.geo.canvas.delete(self.force_vector)
            a = self.rope_points_canvas[-1]
            b = self.tension_direction/scale
            self.force_vector = self.geo.canvas.create_line(a[0],
                                                            a[1],
                                                            a[0] + b[0],
                                                            a[1] - b[1],
                                                            fill='green',arrow=tkr.LAST,width=3)


    def get_stretch(self):
        length = self.compute_length()
        self.stretch = max(0, length - self.length)

    def compute_tension(self):
        direction = self.compute_tension_direction()
        v         = self.climber.get_speed()
        y         = self.stretch * 2.6/self.length
        vt        = np.dot(v, direction)*(y>0)
        self.tension   = max(0,self.k1 * y + self.k2 * y**3 - self.c * vt * (vt>0))
        self.tension_direction = np.array(self.tension * direction)


if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser('climb')
    parser.add_argument('--dt', type=float, default=0.01, required=False)
    parser.add_argument('--tmax', type=float, default=5, required=False)
    parser.add_argument('--offset_time', action='store_true')
    parser.add_argument('--mode', type=str, default='standard', required=False, choices=['standard','fall','crazy'])
    parser.add_argument('--mass', type=float, default=80, required=False)
    parser.add_argument('--animation', action='store_true')

    args = parser.parse_args()

    g=Geometry(W=10,H=20,scale=2)
    a_grav = 9.8

    mode = args.mode
    mass = args.mass

    if  mode=='standard':
        # standard Tests:
        x           = 0.0
        y0 = y      = 2.3
        anchor_x    = -0.3
        anchor_y    = 0.0
        quickdraw_x = 0.0
        quickdraw_y = 0.0
    elif  mode=='fall':
        # realistic fall
        anchor_x    = -0.3
        anchor_y    = -10.0
        quickdraw_x = 0.0
        quickdraw_y = 5.0
        x           = quickdraw_x + 1.0
        y0 = y      = quickdraw_y + 2.0
    elif mode == 'crazy':
        # crazy whip
        anchor_x = -0.3
        anchor_y = -10.0
        quickdraw_x = 0.0
        quickdraw_y = 5.0
        x = quickdraw_x + 1.0
        y0 = y = quickdraw_y + 5.0
        pass

    c=Climber(g,mass,x,y)
    Fgrav = np.array([0.0,-mass*a_grav])
    c.set_force(Fgrav)

    r=Rope(g)
    r.add_anchor(anchor_x,anchor_y,angle = -90)
    r.add_quickdraw(quickdraw_x,quickdraw_y,-45)
    r.add_climber(c)
    length = r.compute_length()
    r.set_rope_length(length)

    dt     = args.dt
    tstart = 0
    t      = 0
    Time   = []
    Y      = []
    VY     = []
    FX     = []
    FY     = []
    T      = []
    S      = []

    #r.set_constants(self, k1, k2, c)

    r.update()

    p    = c.get_position()
    v    = c.get_speed()
    f    = c.get_force()

    if args.animation:
        g.DrawReferenceGrid(10)
        r.draw()
    start_fall=0.0
    set_tstart=False
    while y > g.bottom and tstart<args.tmax:
        if (args.offset_time and r.tension>0) or not args.offset_time:
            Time.append(tstart)
            Y.append(p[1])
            VY.append(v[1])
            FX.append(f[0])
            FY.append(f[1])
            T.append(r.tension)

            if r.tension>0 and not set_tstart:
                start_fall = tstart
                set_tstart = True

            S.append(r.stretch)
            tstart+=dt

        t+=dt
        c.update(dt)
        r.update()

        F = r.tension_direction + Fgrav
        c.set_force(F)

        p  = c.get_position()
        v = c.get_speed()
        f = c.get_force()

        if args.animation:
            r.draw(scale=50)
            c.draw_force(scale=50)
            g.DisplayTime('%.3f' % t)
            g.tk.update()
            time.sleep(0.01)
        pass

    Time  = np.array(Time)
    Y     = np.array(Y)
    VY    = np.array(VY)
    FX    = np.array(FX)
    FY = np.array(FY)
    T = np.array(T)

    S    = np.array(S)

    Tmax  = T.max()
    Fmax  = FY.max()
    Smax  = S.max()

    print ('*---------------------*')
    print ('Maximum Force.......: %.1f kN' % (Fmax/1000.))
    print ('Maximum Tension.....: %.1f kN' % (Tmax/1000.))
    print ('Maximum Elongation..: %.1f %%' % (Smax/length*100.))
    print ('Maximum Stretch.....: %.1f m' % (Smax))
    print ('Length of rope......: %.1f m' % (length))
    print ('Fall factor.........: %.1f m' % (2*y0/length))
    print ('*---------------------*')

    fig=plt.figure(figsize=(7,21))
    fig.add_subplot(3,1,1)
    plt.plot(Time,Y, label = 'rope simulation')
    plt.plot(Time,y0 - 0.5 * 9.8 * Time**2 , label = 'free falling')
    plt.ylabel('Y[m]')
    plt.legend()

    fig.add_subplot(3,1,2)
    plt.plot(Time,FX,':',label='climber-Fx')
    plt.plot(Time,FY,':',label='climber-Fy')
    plt.plot(Time,T,label='rope tension')

    data = RopeData()
    plt.plot(data[0]+start_fall, data[1], 'gray',label='Edelrid Cobra 10.3mm')

    plt.ylabel('Tension F$_{y}$ [N]')
    plt.xlabel('Time [s]')
    plt.legend()

    fig.add_subplot(3, 1, 3)
    plt.plot(-Y, T)
    plt.ylabel('F$_{y}$ [N]')
    plt.xlabel('Y[m]')

    plt.show()
    g.mainloop()
