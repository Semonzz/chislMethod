import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

def system(x, y):
    dx = 4.0 - 4.0 * x - 2.0 * y
    dy = x * y
    return dx, dy

def analyze_fixed_point(x0, y0):
    a = -4.0   # dP/dx
    b = -2.0   # dP/dy
    c = y0     # dQ/dx
    d = x0     # dQ/dy
    trace = a + d
    det = a * d - b * c
    disc = trace * trace - 4.0 * det
    
    EPS = 1e-12
    if disc > EPS:
        l1 = (trace + np.sqrt(disc)) / 2.0
        l2 = (trace - np.sqrt(disc)) / 2.0
        if l1 * l2 < 0:
            return "седло (неустойчиво)"
        elif l1 > 0 and l2 > 0:
            return "неустойчивый узел"
        elif l1 < 0 and l2 < 0:
            return "устойчивый узел"
        else:
            return "узел (граничный случай)"
    elif abs(disc) < EPS:
        l = trace / 2.0
        if l > 0:
            return "неустойчивый вырожденный узел"
        elif l < 0:
            return "устойчивый вырожденный узел"
        else:
            return "центр / седло-узел?"
    else:
        real = trace / 2.0
        if abs(real) < EPS:
            return "центр"
        elif real > 0:
            return "неустойчивый фокус"
        else:
            return "устойчивый фокус"

def rk3_step(x, y, h):
    dx1, dy1 = system(x, y)
    q1x = h * dx1
    q1y = h * dy1
    
    dx2, dy2 = system(x + q1x / 2.0, y + q1y / 2.0)
    q2x = h * dx2
    q2y = h * dy2
    
    dx3, dy3 = system(x - q1x + 2.0 * q2x, y - q1y + 2.0 * q2y)
    q3x = h * dx3
    q3y = h * dy3
    
    x_new = x + (q1x + 4.0 * q2x + q3x) / 6.0
    y_new = y + (q1y + 4.0 * q2y + q3y) / 6.0
    return x_new, y_new

def integrate_trajectory(x0, y0, t_max, h):
    t = 0.0
    x, y = x0, y0
    xs = [x]
    ys = [y]
    
    while t < t_max:
        x, y = rk3_step(x, y, h)
        t += h
        if abs(x) > 100 or abs(y) > 100:
            break
        xs.append(x)
        ys.append(y)
    
    print(f"Траектория из ({x0},{y0}) закончилась в ({x:.4f}, {y:.4f})")
    return xs, ys

def plot_phase_portrait(x0, y0, point_type):
    size_shape = 0.1
    eps = 0.1
    
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='white', linewidth=0.5)
    ax.axvline(x=0, color='white', linewidth=0.5)
    
    colors = ['red', 'blue', 'green', 'orange']
    initials = [
        [x0 + eps, y0], [x0 - eps, y0],
        [x0, y0 + eps], [x0, y0 - eps]
    ]
    
    for i in range(4):
        xs, ys = integrate_trajectory(initials[i][0], initials[i][1], 6.0, 0.02)
        
        ax.plot(xs, ys, color=colors[i], linewidth=1.5, label=f'Траектория {i+1}')
        
        circle = Circle((xs[0], ys[0]), size_shape * 0.1, color=colors[i])
        ax.add_patch(circle)
        
        square = Rectangle((xs[-1] - size_shape * 0.1, ys[-1] - size_shape * 0.1), 
                          size_shape * 0.2, size_shape * 0.2, color=colors[i])
        ax.add_patch(square)
    
    fixed_point = Circle((x0, y0), size_shape * 0.1, color='yellow', zorder=5)
    ax.add_patch(fixed_point)
    
    ax.set_title(f'Фазовый портрет: особая точка ({x0}, {y0})\nТип: {point_type}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'phase_portrait_{x0}_{y0}.png', dpi=150, bbox_inches='tight')
    plt.show()

def main():
    print("Лабораторная работа №5. Вариант 6")
    print("Система:")
    print("  dx/dt = 4 - 4x - 2y")
    print("  dy/dt = x*y")
    print("Численный метод: Рунге-Кутта 3-го порядка, h = 0.02\n")
    
    points = [(0.0, 2.0), (1.0, 0.0)]
    
    print("Особые точки и их тип (аналитически):")
    for p in points:
        typ = analyze_fixed_point(p[0], p[1])
        print(f"  ({p[0]}, {p[1]}) : {typ}")
    print()
    
    for p in points:
        print(f"\nОсобая точка ({p[0]}, {p[1]})")
        typ = analyze_fixed_point(p[0], p[1])
        print(f"Ожидаемый тип: {typ}")
        plot_phase_portrait(p[0], p[1], typ)
        print("  Сгенерировано 4 траектории для этой точки.")
    
    print("\nЧисленное интегрирование завершено. Все графики сохранены в PNG.")
    print("Аналитические выводы и численные результаты должны совпадать.")

if __name__ == "__main__":
    main()
