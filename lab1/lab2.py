def print_vector(a, name):
    print(name)

    for i in range(len(a)):
        print(f"{name}{i + 1} = {a[i]:.6f}")

    print()


def check_stability(a, b, c):
    n = len(b)

    for i in range(1, n - 1):
        if a[i] == 0 or c[i] == 0:
            return False

    strict = False

    for i in range(n):
        if abs(b[i]) < abs(a[i]) + abs(c[i]):
            return False

        if abs(b[i]) > abs(a[i]) + abs(c[i]):
            strict = True

    return strict


def thomas_algorithm(a, b, c, d):
    n = len(b)

    P = [0.0] * n
    Q = [0.0] * n

    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    for i in range(1, n - 1):
        denominator = b[i] + a[i] * P[i - 1]

        P[i] = -c[i] / denominator
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denominator

    P[n - 1] = 0.0

    denominator = b[n - 1] + a[n - 1] * P[n - 2]
    Q[n - 1] = (d[n - 1] - a[n - 1] * Q[n - 2]) / denominator

    x = [0.0] * n

    x[n - 1] = Q[n - 1]

    for i in range(n - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]

    return P, Q, x


def main():
    with open("2.txt", "r", encoding="utf-8") as file:
        data = file.read().split()

    p = 0

    n = int(data[p])
    p += 1

    a = []

    for i in range(n):
        a.append(float(data[p]))
        p += 1

    b = []

    for i in range(n):
        b.append(float(data[p]))
        p += 1

    c = []

    for i in range(n):
        c.append(float(data[p]))
        p += 1

    d = []

    for i in range(n):
        d.append(float(data[p]))
        p += 1

    if check_stability(a, b, c):
        print("Условие устойчивости и диагонального преобладания выполнено.")
        print()

    else:
        print("Условие устойчивости и диагонального преобладания не выполнено.")

    P, Q, x = thomas_algorithm(a, b, c, d)

    print_vector(P, "P")
    print_vector(Q, "Q")
    print_vector(x, "x")


main()