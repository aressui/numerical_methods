import matplotlib.pyplot as plt


def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    if len(lines) != 3:
        raise ValueError("В файле должны быть три строки")

    x_star = float(lines[0])
    xs = [float(value) for value in lines[1].split()]
    ys = [float(value) for value in lines[2].split()]

    return x_star, xs, ys


def solve_tridiagonal(lower, diag, upper, rhs):
    m = len(diag)

    alpha = [0.0] * m
    beta = [0.0] * m

    if abs(diag[0]) < 1e-12:
        raise ValueError("Нулевой элемент на диагонали")

    alpha[0] = -upper[0] / diag[0]
    beta[0] = rhs[0] / diag[0]

    for i in range(1, m):
        denominator = diag[i] + lower[i] * alpha[i - 1]

        if abs(denominator) < 1e-12:
            raise ValueError("Прогонка неустойчива: нулевой знаменатель")

        alpha[i] = -upper[i] / denominator
        beta[i] = (rhs[i] - lower[i] * beta[i - 1]) / denominator

    x = [0.0] * m
    x[m - 1] = beta[m - 1]

    for i in range(m - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]

    return x


def build_spline(xs, ys):
    n = len(xs) - 1

    if n < 2:
        raise ValueError("Для сплайна нужно минимум 3 узла")

    h = [xs[j + 1] - xs[j] for j in range(n)]

    for step in h:
        if step <= 0:
            raise ValueError("Узлы должны быть упорядочены по возрастанию")

    m = n - 1
    lower = [0.0] * m
    diag = [0.0] * m
    upper = [0.0] * m
    rhs = [0.0] * m

    for j in range(1, n):
        r = j - 1
        lower[r] = h[j - 1]
        diag[r] = 2.0 * (h[j - 1] + h[j])
        upper[r] = h[j]
        rhs[r] = 3.0 * ((ys[j + 1] - ys[j]) / h[j] - (ys[j] - ys[j - 1]) / h[j - 1])

    lower[0] = 0.0
    upper[m - 1] = 0.0

    inner = solve_tridiagonal(lower, diag, upper, rhs)
    c_all = [0.0] + inner + [0.0]

    a = []
    b = []
    c = []
    d = []

    for j in range(n):
        a.append(ys[j])
        b.append((ys[j + 1] - ys[j]) / h[j] - h[j] * (c_all[j + 1] + 2.0 * c_all[j]) / 3.0)
        c.append(c_all[j])
        d.append((c_all[j + 1] - c_all[j]) / (3.0 * h[j]))

    return h, a, b, c, d


def find_interval(xs, x):
    for j in range(len(xs) - 1):
        if xs[j] <= x <= xs[j + 1]:
            return j

    raise ValueError("Точка x* лежит вне отрезка [{}; {}]".format(xs[0], xs[-1]))


def spline_value(xs, a, b, c, d, j, x):
    t = x - xs[j]
    return a[j] + b[j] * t + c[j] * t ** 2 + d[j] * t ** 3


def second_derivative(xs, c, d, j, x):
    t = x - xs[j]
    return 2.0 * c[j] + 6.0 * d[j] * t


def plot_spline(xs, ys, a, b, c, d, x_star, s_star):
    plt.figure(figsize=(9, 5))

    for j in range(len(xs) - 1):
        count = 100
        x_values = [xs[j] + (xs[j + 1] - xs[j]) * i / count for i in range(count + 1)]
        y_values = [spline_value(xs, a, b, c, d, j, x) for x in x_values]
        label = "Участок {}: [{:.1f}; {:.1f}]".format(j + 1, xs[j], xs[j + 1])

        plt.plot(x_values, y_values, linewidth=2, label=label)

    point_label = "(x*, S(x*)) = ({}; {:.5f})".format(x_star, s_star)

    plt.scatter(xs, ys, color="black", zorder=5, label="Узлы")
    plt.scatter([x_star], [s_star], color="red", marker="o", s=60, zorder=6, label=point_label)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Кубический сплайн")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    try:
        x_star, xs, ys = read_input("2.txt")
    except (OSError, ValueError, IndexError):
        print("Ошибка чтения файла.")
        return

    if len(xs) != len(ys):
        print("Ошибка: число узлов и значений функции не совпадает.")
        return

    if len(xs) < 3:
        print("Ошибка: для кубического сплайна нужно минимум три узла.")
        return

    try:
        h, a, b, c, d = build_spline(xs, ys)
        j_star = find_interval(xs, x_star)
    except ValueError as error:
        print("Ошибка:", error)
        return

    n = len(xs) - 1

    print("Исходные узлы и значения функции:")
    print("i   x_i      f_i")
    for i in range(n + 1):
        print("{}   {:.4f}   {:.5f}".format(i, xs[i], ys[i]))

    print("\nКоэффициенты кубического сплайна:")
    print("S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3")
    print("i   интервал       a_i          b_i          c_i          d_i")
    for j in range(n):
        print("{}   [{:.1f}; {:.1f}]  {:10.6f}   {:10.6f}   {:10.6f}   {:10.6f}".format(
            j + 1, xs[j], xs[j + 1], a[j], b[j], c[j], d[j]))

    s_star = spline_value(xs, a, b, c, d, j_star, x_star)

    print("\nТочка x* = {}".format(x_star))
    print("Номер интервала, содержащего x*: {} ([{:.4f}; {:.4f}])".format(j_star + 1, xs[j_star], xs[j_star + 1]))
    print("S(x*) = {:.6f}".format(s_star))

    print("\nПроверка прохождения сплайна через узлы:")
    print("i   x_i      f_i        S(x_i)      |S(x_i) - f_i|")
    for i in range(n + 1):
        j = i - 1
        if i == 0:
            j = 0

        value = spline_value(xs, a, b, c, d, j, xs[i])
        print("{}   {:.4f}   {:.5f}   {:.6f}   {:.2e}".format(i, xs[i], ys[i], value, abs(value - ys[i])))

    print("\nВторая производная на концах (граничные условия):")
    print("S''(x_0) = {:.6f}".format(second_derivative(xs, c, d, 0, xs[0])))
    print("S''(x_n) = {:.6f}".format(second_derivative(xs, c, d, n - 1, xs[n])))

    plot_spline(xs, ys, a, b, c, d, x_star, s_star)


if __name__ == "__main__":
    main()