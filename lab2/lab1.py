
import math
import matplotlib.pyplot as plt


EPSILON = 0.0001
MAX_ITERATIONS = 10000


def f(x):
    return math.cos(x) + 0.25 * x - 0.5


def df(x):
    return -math.sin(x) + 0.25


def d2f(x):
    return -math.cos(x)


def phi(x):
    return math.acos(0.5 - 0.25 * x)


def dphi(x):
    t = 0.5 - 0.25 * x
    return 0.25 / math.sqrt(1.0 - t * t)


def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    scan = [float(value) for value in lines[0].split()]
    segment = [float(value) for value in lines[1].split()]
    x0_newton = float(lines[2])
    x0_iter = float(lines[3])

    return scan, segment, x0_newton, x0_iter


def localize(start, end, step):
    print("Табулирование функции:")
    print("{:>10} {:>14}".format("x", "f(x)"))

    segments = []
    n = int(round((end - start) / step))
    prev_x = start
    prev_f = f(start)

    print("{:10.4f} {:14.6f}".format(prev_x, prev_f))

    for i in range(1, n + 1):
        x = start + i * step
        value = f(x)

        print("{:10.4f} {:14.6f}".format(x, value))

        if prev_f * value < 0:
            segments.append((prev_x, x))

        prev_x = x
        prev_f = value

    print("\nОтрезки со сменой знака:")
    for left, right in segments:
        print("[{:.4f}; {:.4f}]".format(left, right))

    return segments


def newton(x0, epsilon):
    x = x0

    print("{:>4} {:>14} {:>14} {:>14} {:>14}".format(
        "k", "x_k", "f(x_k)", "f'(x_k)", "delta x"
    ))

    for k in range(MAX_ITERATIONS):
        fx = f(x)
        derivative = df(x)

        if abs(derivative) < 1e-12:
            raise ValueError("Производная близка к нулю")

        x_new = x - fx / derivative
        diff = abs(x_new - x)

        print("{:4d} {:14.8f} {:14.8f} {:14.8f} {:14.8f}".format(
            k, x, fx, derivative, diff
        ))

        if diff < epsilon:
            return x_new, k + 1

        x = x_new

    raise ValueError("Метод Ньютона не сошелся")


def simple_iteration(x0, a, b, epsilon):
    samples = 1000
    q = 0.0
    inside = True

    for i in range(samples + 1):
        x = a + (b - a) * i / samples
        q = max(q, abs(dphi(x)))

        if not a <= phi(x) <= b:
            inside = False

    print("\nПроверка условий сходимости:")
    print("phi([a, b]) принадлежит [a, b] —", "да" if inside else "нет")
    print("Коэффициент сходимости q = {:.6f}".format(q))

    if not inside or q >= 1:
        print("Условия сходимости не выполнены")
        return x0, 0

    x = x0

    print("{:>4} {:>14} {:>14} {:>14}".format(
        "k", "x_k", "x_(k+1)", "delta x"
    ))

    for k in range(MAX_ITERATIONS):
        x_new = phi(x)
        diff = abs(x_new - x)

        print("{:4d} {:14.8f} {:14.8f} {:14.8f}".format(
            k, x, x_new, diff
        ))

        if q / (1.0 - q) * diff <= epsilon and abs(f(x_new)) <= epsilon:
            return x_new, k + 1

        x = x_new

    raise ValueError("Метод простой итерации не сошелся")


def plot_function(a, b):
    count = 1000
    x_values = [-2 + i * 0.01 for i in range(count + 1)]

    y_cos = [math.cos(x) for x in x_values]
    y_line = [0.5 - 0.25 * x for x in x_values]

    plt.figure(figsize=(9, 5))
    plt.plot(x_values, y_cos, label="y = cos(x)", linewidth=2)
    plt.plot(x_values, y_line, label="y = 0.5 - 0.25x", linewidth=2)
    plt.axvspan(a, b, alpha=0.15, label="Отрезок поиска")
    plt.xlim(-2, 8)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Графическое решение уравнения")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    scan, segment, x0_newton, x0_iter = read_input("1.txt")
    a, b = segment

    localize(scan[0], scan[1], scan[2])

    print("\nОтрезок поиска: [{:.4f}; {:.4f}]".format(a, b))
    print("f(a) * f(b) = {:.8f}".format(f(a) * f(b)))
    print("f(x0) * f''(x0) = {:.8f}".format(
        f(x0_newton) * d2f(x0_newton)
    ))

    if f(a) * f(b) >= 0:
        print("На концах отрезка нет смены знака.")
        return

    print("\nМетод Ньютона")
    try:
        root_n, k_n = newton(x0_newton, EPSILON)
        print("\nКорень: {:.8f}".format(root_n))
        print("Значение функции в найденном корне: {:.8f}".format(f(root_n)))
        print("Итераций:", k_n)
    except ValueError as error:
        print("Ошибка:", error)

    print("\nМетод простой итерации")
    root_i, k_i = simple_iteration(x0_iter, a, b, EPSILON)

    if k_i > 0:
        print("\nКорень: {:.8f}".format(root_i))
        print("Значение функции в найденном корне: {:.8f}".format(f(root_i)))
        print("Итераций:", k_i)

    plot_function(a, b)


if __name__ == "__main__":
    main()
