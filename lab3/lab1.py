import math
import matplotlib.pyplot as plt


def f(x):
    return math.log(x)


def d4f(x):
    return -6.0 / x ** 4


def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    if len(lines) != 3:
        raise ValueError("В файле должны быть три строки")

    x_star = float(lines[0])
    node_sets = [[float(value) for value in line.split()] for line in lines[1:]]

    return x_star, node_sets


def poly_multiply(p, q):
    result = [0.0] * (len(p) + len(q) - 1)

    for i in range(len(p)):
        for j in range(len(q)):
            result[i + j] += p[i] * q[j]

    return result


def format_polynomial(coefficients):
    text = ""

    for power in range(len(coefficients) - 1, -1, -1):
        value = abs(coefficients[power])

        if power == 0:
            term = "{:.6f}".format(value)
        elif power == 1:
            term = "{:.6f}*x".format(value)
        else:
            term = "{:.6f}*x^{}".format(value, power)

        if power == len(coefficients) - 1:
            if coefficients[power] < 0:
                text += "-" + term
            else:
                text += term
        else:
            if coefficients[power] < 0:
                text += " - " + term
            else:
                text += " + " + term

    return text


def lagrange_value(xs, ys, x):
    total = 0.0

    for i in range(len(xs)):
        term = ys[i]

        for j in range(len(xs)):
            if j != i:
                term *= (x - xs[j]) / (xs[i] - xs[j])

        total += term

    return total


def lagrange_coefficients(xs, ys):
    n = len(xs)
    result = [0.0] * n

    for i in range(n):
        basis = [1.0]
        denominator = 1.0

        for j in range(n):
            if j != i:
                basis = poly_multiply(basis, [-xs[j], 1.0])
                denominator *= xs[i] - xs[j]

        for k in range(n):
            result[k] += ys[i] * basis[k] / denominator

    return result


def divided_differences(xs, ys):
    n = len(xs)
    table = [ys[:]]

    for k in range(1, n):
        previous = table[k - 1]
        row = []

        for i in range(n - k):
            row.append((previous[i + 1] - previous[i]) / (xs[i + k] - xs[i]))

        table.append(row)

    return table


def print_divided_differences(xs, table):
    n = len(xs)

    header = "{:>8} {:>12}".format("x_i", "f(x_i)")
    for k in range(1, n):
        header += " {:>12}".format(str(k) + " порядок")
    print(header)

    for i in range(n):
        line = "{:8.4f} {:12.6f}".format(xs[i], table[0][i])

        for k in range(1, n - i):
            line += " {:12.6f}".format(table[k][i])

        print(line)


def newton_value(xs, coefficients, x):
    result = coefficients[-1]

    for k in range(len(coefficients) - 2, -1, -1):
        result = result * (x - xs[k]) + coefficients[k]

    return result


def format_newton(xs, coefficients):
    text = "{:.6f}".format(coefficients[0])
    factors = []

    for k in range(1, len(coefficients)):
        factors.append("(x - {:.4f})".format(xs[k - 1]))

        sign = "+"
        if coefficients[k] < 0:
            sign = "-"

        text += " {} {:.6f}*{}".format(sign, abs(coefficients[k]), "*".join(factors))

    return text


def newton_expanded(xs, coefficients):
    n = len(xs)
    result = [0.0] * n
    basis = [1.0]

    for k in range(n):
        for i in range(len(basis)):
            result[i] += coefficients[k] * basis[i]

        basis = poly_multiply(basis, [-xs[k], 1.0])

    return result


def error_bound(xs, x):
    n = len(xs)
    left = min(xs)
    right = max(xs)
    samples = 1000

    m = 0.0
    for i in range(samples + 1):
        point = left + (right - left) * i / samples
        m = max(m, abs(d4f(point)))

    omega = 1.0
    for node in xs:
        omega *= x - node

    return m, m / math.factorial(n) * abs(omega)


def plot_result(xs, ys, coefficients, x_star, title):
    left = max(min(xs) - 0.1, 0.05)
    right = max(xs) + 0.1
    count = 500

    x_values = [left + (right - left) * i / count for i in range(count + 1)]
    y_exact = [f(x) for x in x_values]
    y_lagrange = [lagrange_value(xs, ys, x) for x in x_values]
    y_newton = [newton_value(xs, coefficients, x) for x in x_values]

    plt.figure(figsize=(9, 5))
    plt.plot(x_values, y_exact, label="y = ln(x)", linewidth=2)
    plt.plot(x_values, y_lagrange, label="Многочлен Лагранжа", linewidth=2)
    plt.plot(x_values, y_newton, "--", label="Многочлен Ньютона", linewidth=2)
    lagrange_star = lagrange_value(xs, ys, x_star)

    plt.scatter(xs, ys, color="black", zorder=5, label="Узлы интерполяции")
    plt.scatter([x_star], [lagrange_star], color="red", marker="o", s=60,
                zorder=6, label="Точка x* = {:.2f}".format(x_star))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title + ", y = ln(x)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def solve(xs, x_star, title):
    n = len(xs)

    if n < 2:
        raise ValueError("Для интерполяции нужно минимум два узла")

    if x_star <= 0:
        raise ValueError("Точка x* должна быть положительной")

    for i in range(n):
        if xs[i] <= 0:
            raise ValueError("Для ln x все узлы должны быть положительными")

        for j in range(i + 1, n):
            if xs[i] == xs[j]:
                raise ValueError("Узлы интерполяции должны быть различными")

    ys = [f(x) for x in xs]

    print("\n" + title)

    print("\nУзлы интерполяции:")
    print("i   x_i      y_i")
    for i in range(n):
        print("{}   {:.4f}   {:.6f}".format(i, xs[i], ys[i]))

    lagrange_poly = lagrange_coefficients(xs, ys)

    print("\nМногочлен Лагранжа:")
    print("L(x) =", format_polynomial(lagrange_poly))

    table = divided_differences(xs, ys)
    newton_coefficients = [table[k][0] for k in range(n)]

    print("\nТаблица разделённых разностей:")
    print_divided_differences(xs, table)

    print("\nМногочлен Ньютона:")
    print("N(x) =", format_newton(xs, newton_coefficients))
    print("\nМногочлен Ньютона (после раскрытия скобок):")
    print("N(x) =", format_polynomial(newton_expanded(xs, newton_coefficients)))

    l_value = lagrange_value(xs, ys, x_star)
    n_value = newton_value(xs, newton_coefficients, x_star)
    exact = f(x_star)

    print("\nРезультаты в точке x* = {}:".format(x_star))
    print("L(x*) = {:.6f}".format(l_value))
    print("N(x*) = {:.6f}".format(n_value))
    print("f(x*) = {:.6f}".format(exact))
    print("|f(x*) - L(x*)| = {:.6f}".format(abs(exact - l_value)))
    print("|f(x*) - N(x*)| = {:.6f}".format(abs(exact - n_value)))
    print("|L(x*) - N(x*)| = {:.2e}".format(abs(l_value - n_value)))

    if n == 4:
        m, bound = error_bound(xs, x_star)
        print("\nОценка погрешности по формуле (3.9): {:.6f}".format(bound))

    plot_result(xs, ys, newton_coefficients, x_star, title)


def main():
    try:
        x_star, node_sets = read_input("1.txt")
    except (OSError, ValueError, IndexError):
        print("Ошибка чтения файла.")
        return

    names = ["а", "б"]

    for index in range(len(node_sets)):
        if index < len(names):
            name = names[index]
        else:
            name = str(index + 1)

        title = "Набор узлов {}:".format(name)

        try:
            solve(node_sets[index], x_star, title)
        except ValueError as error:
            print("Ошибка:", error)


if __name__ == "__main__":
    main()