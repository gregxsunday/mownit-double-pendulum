import math
from math import sin, cos
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import tkinter as tk

N = 4

def derivs(xin, yin, properties):
    L1, L2 = properties['lengths']
    M1, M2 = properties['masses']
    G = properties['g']

    dydx = [0 for _ in range(N)]
    dydx[0] = yin[1]

    de = yin[2] - yin[0]

    den1 = (M1+M2)*L1 - M2*L1*cos(de)**2
    dydx[1] = (M2*L1*yin[1]*yin[1]*sin(de)*cos(de) + M2*G*sin(yin[2])*cos(de) + M2*L2*yin[3]*yin[3]*sin(de) - (M1+M2)*G*sin(yin[0]))/den1

    dydx[2] = yin[3]
    den2 = (L2 / L1) * den1

    dydx[3] = (-M2 * L2 * yin[3] * yin[3] * sin(de) * cos(de) + (M1 + M2) * G * sin(yin[0]) * cos(de) - (M1 + M2) * L1 * yin[1] * yin[1] * sin(de) - (M1 + M2) * G * sin(yin[2]))/den2
    return dydx


def runge_kutta(xin, yin, h, properties):
    yout = [0 for _ in range(N)]
    hh = 0.5*h
    xh = xin + h
    yt = [0 for _ in range(N)]
    k = [[0 for _ in range(N)] for _ in range(N)]

    dydx = derivs(xin, yin, properties)
    for i in range(N):
        k[0][i] = h*dydx[i]
        yt[i] = yin[i] + 0.5*k[0][i]

    dydxt = derivs(xh, yt, properties)
    for i in range(N):
        k[1][i] = h * dydxt[i]
        yt[i] = yin[i] + 0.5 * k[1][i]

    dydxt = derivs(xh, yt, properties)
    for i in range(N):
        k[2][i] = h * dydxt[i]
        yt[i] = yin[i] + k[2][i]

    dydxt = derivs(xin + h, yt, properties)
    for i in range(N):
        k[3][i] = h * dydxt[i]
        yout[i] = yin[i] + k[0][i]/6 + k[1][i]/3 + k[2][i]/3 + k[3][i]/6

    return yout


def show_animation(time_start=0, time_end=10, initial_theta1=30, initial_theta2=60, g=9.81, speed=1, len1=1, len2=1, m1=1, m2=1):
    speed = int(speed)
    duration = int(time_end) - int(time_start)
    fps = 30
    initial_thetas = [int(initial_theta1), int(initial_theta2)]
    initial_angular_velocities = [0, 0]
    steps = duration*fps
    lengths = [float(len1), float(len2)]
    masses = [float(m1), float(m2)]

    properties = {'lengths': lengths, 'masses': masses, 'g': float(g)}
    h = duration / (steps - 1)
    t = [int(time_start) + i*h for i in range(steps)]

    thetas = [[0 for _ in range(steps)] for _ in range(2)]
    angular_velocities = [[0 for _ in range(steps)] for _ in range(2)]

    x1, y1, x2, y2 = [], [], [], []

    for i in range(2):
        thetas[i][0] = initial_thetas[i]*math.pi/180
        angular_velocities[i][0] = initial_angular_velocities[i]*math.pi/180

    for i in range(steps - 1):
        x1.append(lengths[0]*sin(thetas[0][i]))
        y1.append(-lengths[0]*cos(thetas[0][i]))
        x2.append(x1[i] + lengths[1]*sin(thetas[1][i]))
        y2.append(y1[i] - lengths[1]*cos(thetas[1][i]))

        yin = [thetas[0][i], angular_velocities[0][i], thetas[1][i], angular_velocities[1][i]]
        yout = runge_kutta(t[i], yin, h, properties)
        thetas[0][i + 1], angular_velocities[0][i + 1], thetas[1][i + 1], angular_velocities[1][i + 1] = yout

    fig = plt.figure()
    fig.canvas.set_window_title('Double Pendulum')
    plot_size = lengths[0] + lengths[1] + 0.5
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-plot_size, plot_size), ylim=(-plot_size, plot_size))
    ax.set_aspect('equal')
    ax.grid()
    print(type(ax))
    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text


    def animate(i):
        thisx = [0, x1[i], x2[i]]
        thisy = [0, y1[i], y2[i]]

        line.set_data(thisx, thisy)
        time_text.set_text(time_template % (t[i]))
        return line, time_text

    ani = animation.FuncAnimation(fig, animate, range(steps - 1), interval=h*1000*(1/speed), blit=True, init_func=init)
    plt.show()


def init_gui():
    root = tk.Tk()
    root.title('Double pendulum parameters')

    properties = {
        'time start [s]': 0,
        'time end [s]': 10,
        'initial θ\u2081 [\u00b0]': 30,
        'initial θ\u2082 [\u00b0]': 60,
        'g [m/s\u00b2]': 9.81,
        'speed': 1,
        'L\u2081 [m]': 1,
        'L\u2082 [m]': 1,
        'M\u2081 [kg]': 1,
        'M\u2082 [kg]': 1
    }
    e = []
    for i, (key, value) in enumerate(properties.items()):
        tk.Label(root, text=key).grid(row=i)
        e.append(tk.Entry(root))
        e[i].insert(0, str(value))
        e[i].grid(row=i, column=1)

    def callback():
        props = list(map(lambda e: e.get(), e))
        show_animation(*props)

    btn = tk.Button(root, text='Animate', command=callback)
    btn.grid(row=len(properties), column=1)
    root.mainloop()


if __name__ == '__main__':
    init_gui()