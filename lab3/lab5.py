import math


def f(x):
    return 1.0 / ((2 * x + 7) * (3 * x + 4))


def exact_integral(a, b):
    upper = math.log(abs(3 * b + 4)) - math.log(abs(2 * b + 7))
    lower = math.log(abs(3 * a + 4)) - math.log(abs(2 * a + 7))
    return (upper - lower) / 13


def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    a, b = [float(value) for value in lines[0].split()]
    h1, h2 = [float(value) for value in lines[1].split()]

    return a, b, h1, h2


def steps_count(a, b, h):
    n = round((b - a) / h)

    if n < 1 or abs((b - a) / h - n) > 1e-9:
        raise ValueError("Отрезок не делится на целое число шагов h = {}".format(h))

    return n


def left_rectangles(a, b, h):
    n = steps_count(a, b, h)
    total = 0.0

    for i in range(n):
        total += f(a + i * h)

    return h * total


def right_rectangles(a, b, h):
    n = steps_count(a, b, h)
    total = 0.0

    for i in range(1, n + 1):
        total += f(a + i * h)

    return h * total


def middle_rectangles(a, b, h):
    n = steps_count(a, b, h)
    total = 0.0

    for i in range(n):
        total += f(a + (i + 0.5) * h)

    return h * total


def trapezoids(a, b, h):
    n = steps_count(a, b, h)
    total = (f(a) + f(b)) / 2

    for i in range(1, n):
        total += f(a + i * h)

    return h * total


def simpson(a, b, h):
    n = steps_count(a, b, h)

    if n % 2 != 0:
        raise ValueError("Для метода Симпсона число шагов должно быть чётным")

    total = f(a) + f(b)

    for i in range(1, n):
        if i % 2 == 1:
            total += 4 * f(a + i * h)
        else:
            total += 2 * f(a + i * h)

    return h / 3 * total


def main():
    try:
        a, b, h1, h2 = read_input("5.txt")
    except (OSError, ValueError, IndexError):
        print("Ошибка чтения файла.")
        return

    methods = [
        ("Левые прямоугольники", left_rectangles, 1),
        ("Правые прямоугольники", right_rectangles, 1),
        ("Средние прямоугольники", middle_rectangles, 2),
        ("Трапеции", trapezoids, 2),
        ("Симпсон", simpson, 4),
    ]

    k = h1 / h2

    print("Интеграл от a = {} до b = {}, h1 = {}, h2 = {}".format(a, b, h1, h2))
    print()
    print("{:<24} {:>10} {:>10} {:>3} {:>12} {:>12}".format(
        "Метод", "I(h1)", "I(h2)", "p", "Рунге-Ромберг", "Погрешность"))

    results = []

    for name, method, p in methods:
        try:
            value1 = method(a, b, h1)
            value2 = method(a, b, h2)
        except ValueError as error:
            print("{}: {}".format(name, error))
            return

        refined = value2 + (value2 - value1) / (k ** p - 1)
        estimate = abs(value2 - value1) / (k ** p - 1)

        results.append((name, value2, refined))

        print("{:<24} {:10.6f} {:10.6f} {:>3} {:12.6f} {:12.6f}".format(
            name, value1, value2, p, refined, estimate))

    exact = exact_integral(a, b)

    print("\nТочное значение интеграла: {:.6f}".format(exact))
    print()
    print("{:<24} {:>14} {:>14}".format("Метод", "Ошибка I(h2)", "Ошибка уточн."))

    best_plain = 0
    best_refined = 0

    for i in range(len(results)):
        name, value2, refined = results[i]
        error_plain = abs(value2 - exact)
        error_refined = abs(refined - exact)

        print("{:<24} {:14.6f} {:14.6f}".format(name, error_plain, error_refined))

        if error_plain < abs(results[best_plain][1] - exact):
            best_plain = i

        if error_refined < abs(results[best_refined][2] - exact):
            best_refined = i

    print("\nНаиболее точный результат с шагом h2: {}".format(results[best_plain][0]))
    print("Наиболее точный результат после уточнения: {}".format(results[best_refined][0]))


if __name__ == "__main__":
    main()