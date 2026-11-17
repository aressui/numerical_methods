import matplotlib.pyplot as plt


def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    xs = [float(value) for value in lines[0].split()]
    ys = [float(value) for value in lines[1].split()]

    return xs, ys


def solve_gauss(matrix, rhs):
    n = len(rhs)
    m = [matrix[i][:] + [rhs[i]] for i in range(n)]

    for k in range(n):
        p = k
        for i in range(k + 1, n):
            if abs(m[i][k]) > abs(m[p][k]):
                p = i

        if abs(m[p][k]) < 1e-12:
            raise ValueError("Система вырождена")

        m[k], m[p] = m[p], m[k]

        for i in range(k + 1, n):
            factor = m[i][k] / m[k][k]

            for j in range(k, n + 1):
                m[i][j] -= factor * m[k][j]

    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        s = m[i][n]

        for j in range(i + 1, n):
            s -= m[i][j] * x[j]

        x[i] = s / m[i][i]

    return x


def least_squares(xs, ys, degree):
    size = degree + 1
    matrix = []
    rhs = []

    for k in range(size):
        row = []

        for i in range(size):
            s = 0.0
            for x in xs:
                s += x ** (k + i)
            row.append(s)

        matrix.append(row)

        s = 0.0
        for j in range(len(xs)):
            s += ys[j] * xs[j] ** k
        rhs.append(s)

    return solve_gauss(matrix, rhs)


def poly_value(coefficients, x):
    result = 0.0

    for c in reversed(coefficients):
        result = result * x + c

    return result


def mse(xs, ys, coefficients):
    total = 0.0

    for i in range(len(xs)):
        total += (poly_value(coefficients, xs[i]) - ys[i]) ** 2

    return total / len(xs)


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


def plot_result(xs, ys, polynomials):
    left = min(xs) - 0.2
    right = max(xs) + 0.2
    count = 500

    x_values = [left + (right - left) * i / count for i in range(count + 1)]

    plt.figure(figsize=(9, 5))
    plt.scatter(xs, ys, color="black", zorder=5, label="Исходные данные")

    for degree in range(1, len(polynomials) + 1):
        y_values = [poly_value(polynomials[degree - 1], x) for x in x_values]
        plt.plot(x_values, y_values, linewidth=2, label="P{}(x)".format(degree))

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Метод наименьших квадратов")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    try:
        xs, ys = read_input("3.txt")
    except (OSError, ValueError, IndexError):
        print("Ошибка чтения файла.")
        return

    if len(xs) != len(ys):
        print("Ошибка: число значений x и y не совпадает.")
        return

    if len(xs) < 4:
        print("Ошибка: для многочлена третьей степени нужно минимум 4 точки.")
        return

    print("Исходные данные:")
    print("i   x_i      y_i")
    for i in range(len(xs)):
        print("{}   {:.4f}   {:.5f}".format(i, xs[i], ys[i]))

    polynomials = []
    errors = []

    for degree in range(1, 4):
        try:
            coefficients = least_squares(xs, ys, degree)
        except ValueError as error:
            print("Ошибка:", error)
            return

        polynomials.append(coefficients)
        errors.append(mse(xs, ys, coefficients))

        print("\nМногочлен степени {}:".format(degree))

        for i in range(len(coefficients)):
            print("a{} = {:.6f}".format(i, coefficients[i]))

        print("P{}(x) = {}".format(degree, format_polynomial(coefficients)))

    print()
    for i in range(3):
        print("MSE{} = {:.6f}".format(i + 1, errors[i]))

    best = 0
    for i in range(1, 3):
        if errors[i] < errors[best]:
            best = i

    print("\nЛучше всего приближает данные многочлен P{}(x), у него наименьшая ошибка.".format(best + 1))

    plot_result(xs, ys, polynomials)


if __name__ == "__main__":
    main()