def read_input(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    x_star = float(lines[0])
    xs = [float(value) for value in lines[1].split()]
    ys = [float(value) for value in lines[2].split()]

    return x_star, xs, ys


def choose_nodes(xs, x):
    n = len(xs)

    for k in range(n):
        if abs(xs[k] - x) < 1e-9:
            if k == 0:
                return 0, k
            if k == n - 1:
                return n - 3, k
            return k - 1, k

    for i in range(n - 1):
        if xs[i] < x < xs[i + 1]:
            if i + 2 <= n - 1:
                return i, None
            return i - 1, None

    raise ValueError("Точка x* лежит вне таблицы")


def main():
    try:
        x_star, xs, ys = read_input("4.txt")
    except (OSError, ValueError, IndexError):
        print("Ошибка чтения файла.")
        return

    if len(xs) != len(ys) or len(xs) < 3:
        print("Ошибка: в таблице должно быть не меньше 3 узлов.")
        return

    try:
        start, node = choose_nodes(xs, x_star)
    except ValueError as error:
        print("Ошибка:", error)
        return

    print("Исходная таблица:")
    print("i   x_i      y_i")
    for i in range(len(xs)):
        print("{}   {:.4f}   {:.4f}".format(i, xs[i], ys[i]))

    print("\nТочка x* = {}".format(x_star))

    x0, x1, x2 = xs[start], xs[start + 1], xs[start + 2]
    y0, y1, y2 = ys[start], ys[start + 1], ys[start + 2]

    print("Узлы, использованные при вычислении:")
    print("x = {}, {}, {}".format(x0, x1, x2))

    f01 = (y1 - y0) / (x1 - x0)
    f12 = (y2 - y1) / (x2 - x1)
    f012 = (f12 - f01) / (x2 - x0)

    first = f01 + f012 * (2 * x_star - x0 - x1)
    second = 2 * f012

    print("\nПервая производная:")

    if node is not None:
        if node > 0:
            left = (ys[node] - ys[node - 1]) / (xs[node] - xs[node - 1])
            print("левосторонняя: {:.6f}".format(left))

        if node < len(xs) - 1:
            right = (ys[node + 1] - ys[node]) / (xs[node + 1] - xs[node])
            print("правосторонняя: {:.6f}".format(right))

    print("по формуле (3.24): {:.6f}".format(first))

    print("\nВторая производная: {:.6f}".format(second))


if __name__ == "__main__":
    main()