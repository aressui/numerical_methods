import math
import matplotlib.pyplot as plt


MAX_ITERATIONS = 10000


def f1(x1, x2):
    return x1 - math.cos(x2) - 1


def f2(x1, x2):
    return x2 - math.log10(x1 + 1) - 2


def df1_dx1(x1, x2):
    return 1.0


def df1_dx2(x1, x2):
    return math.sin(x2)


def df2_dx1(x1, x2):
    return -1.0 / ((x1 + 1) * math.log(10))


def df2_dx2(x1, x2):
    return 1.0


def phi1(x1, x2):
    return 1 + math.cos(x2)


def phi2(x1, x2):
    return 2 + math.log10(x1 + 1)


def plot_function():
    count = 1000
    x1_values = [2.0 * i / count for i in range(count + 1)]

    x2_curve1 = [math.acos(x1 - 1) for x1 in x1_values]
    x2_curve2 = [2 + math.log10(x1 + 1) for x1 in x1_values]

    plt.figure(figsize=(9, 5))
    plt.plot(x1_values, x2_curve1, label="x1 = 1 + cos(x2)")
    plt.plot(x1_values, x2_curve2, label="x2 = 2 + lg(x1 + 1)")
    plt.scatter([0.5], [2.2], label="Начальное приближение")
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.title("Графическое решение системы")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def newton(x1, x2, epsilon):
    print("\nМетод Ньютона")
    print("{:>2} {:>10} {:>10} {:>11} {:>11} {:>9} {:>9} "
          "{:>10} {:>10} {:>10} {:>10} {:>10}".format(
              "k", "x1", "x2", "f1", "f2", "j11", "j12",
              "j21", "j22", "det(J)", "dx1", "dx2"
          ))

    for k in range(MAX_ITERATIONS):
        value1 = f1(x1, x2)
        value2 = f2(x1, x2)

        j11 = df1_dx1(x1, x2)
        j12 = df1_dx2(x1, x2)
        j21 = df2_dx1(x1, x2)
        j22 = df2_dx2(x1, x2)

        det = j11 * j22 - j12 * j21

        if abs(det) < 1e-12:
            raise ValueError("Определитель матрицы Якоби близок к нулю")

        dx1 = (-value1 * j22 + j12 * value2) / det
        dx2 = (j21 * value1 - j11 * value2) / det

        print("{:2d} {:10.6f} {:10.6f} {:11.3e} {:11.3e} "
              "{:9.5f} {:9.5f} {:10.5f} {:10.5f} "
              "{:10.5f} {:10.3e} {:10.3e}".format(
                  k, x1, x2, value1, value2, j11, j12,
                  j21, j22, det, dx1, dx2
              ))

        x1_new = x1 + dx1
        x2_new = x2 + dx2

        diff = max(abs(x1_new - x1), abs(x2_new - x2))

        if diff <= epsilon:
            return x1_new, x2_new, k + 1

        x1 = x1_new
        x2 = x2_new

    raise ValueError("Метод Ньютона не сошелся")


def check_convergence():
    q = max(math.sin(2.0), 1.0 / math.log(10))

    print("\nПроверка условий сходимости:")
    print("Область: 0 <= x1 <= 2, 2 <= x2 <= 3")
    print("Коэффициент сходимости q = {:.6f}".format(q))

    if q < 1:
        print("Условия сходимости выполнены")
        return q

    print("Условия сходимости не выполнены")
    return None


def simple_iteration(x1, x2, epsilon):
    print("\nМетод простой итерации")

    q = check_convergence()

    if q is None:
        raise ValueError("Условия сходимости не выполнены")

    print("{:>3} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
        "k", "x1", "x2", "x1 новый", "x2 новый", "delta x"
    ))

    for k in range(MAX_ITERATIONS):
        x1_new = phi1(x1, x2)
        x2_new = phi2(x1, x2)

        diff = max(abs(x1_new - x1), abs(x2_new - x2))

        print("{:3d} {:12.7f} {:12.7f} {:12.7f} {:12.7f} {:12.5e}".format(
            k, x1, x2, x1_new, x2_new, diff
        ))

        if q / (1.0 - q) * diff <= epsilon:
            return x1_new, x2_new, k + 1

        x1 = x1_new
        x2 = x2_new

    raise ValueError("Метод простой итерации не сошелся")


def print_result(method, x1, x2, iterations):
    print("\n{}".format(method))
    print("x1 = {:.8f}".format(x1))
    print("x2 = {:.8f}".format(x2))
    print("f1(x1, x2) = {:.18f}".format(f1(x1, x2)))
    print("f2(x1, x2) = {:.18f}".format(f2(x1, x2)))
    print("Количество итераций:", iterations)


def main():
    try:
        with open("2.txt", "r", encoding="utf-8") as file:
            epsilon = float(file.readline())
            x1_0 = float(file.readline())
            x2_0 = float(file.readline())
    except (OSError, ValueError):
        print("Ошибка чтения input.txt. Проверьте файл и его содержимое.")
        return

    if epsilon <= 0 or x1_0 <= -1:
        print("Ошибка: точность должна быть положительной, x1 должен быть больше -1.")
        return

    print("Точность:", epsilon)
    print("Начальное приближение: x1 = {:.2f}, x2 = {:.2f}".format(
        x1_0, x2_0
    ))

    plot_function()

    try:
        root_n1, root_n2, k_n = newton(x1_0, x2_0, epsilon)
        print_result("Результат метода Ньютона", root_n1, root_n2, k_n)
    except ValueError as error:
        print("Ошибка метода Ньютона:", error)
        k_n = None

    try:
        root_i1, root_i2, k_i = simple_iteration(
            x1_0, x2_0, epsilon
        )
        print_result(
            "Результат метода простой итерации",
            root_i1, root_i2, k_i
        )
    except ValueError as error:
        print("Ошибка метода простой итерации:", error)
        k_i = None

    if k_n is not None and k_i is not None:
        print("\nСравнение методов:")
        print("Метод Ньютона: {} итераций".format(k_n))
        print("Метод простой итерации: {} итераций".format(k_i))


if __name__ == "__main__":
    main()