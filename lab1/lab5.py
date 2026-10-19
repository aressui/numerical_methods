import math


def read_input(filename, has_epsilon=True):
    with open(filename, "r") as file:
        lines = file.readlines()

    if has_epsilon:
        epsilon = float(lines[0])
        lines = lines[1:]
    else:
        epsilon = None

    A = []

    for line in lines:
        row = []

        for value in line.split():
            row.append(float(value))

        if len(row) > 0:
            A.append(row)

    return epsilon, A


def copy_matrix(A):
    B = []

    for row in A:
        B.append(row[:])

    return B


def identity_matrix(n):
    I = []

    for i in range(n):
        row = []

        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                row.append(0.0)

        I.append(row)

    return I


def print_matrix(A):
    for row in A:
        for value in row:
            print(round(value, 6), end=" ")

        print()


def multiply(A, B):
    n = len(A)
    C = []

    for i in range(n):
        row = []

        for j in range(n):
            s = 0.0

            for k in range(n):
                s += A[i][k] * B[k][j]

            row.append(s)

        C.append(row)

    return C


def qr_decomposition(A):
    n = len(A)

    R = copy_matrix(A)
    Q = identity_matrix(n)

    for k in range(n - 1):
        norm = 0.0

        for i in range(k, n):
            norm += R[i][k] ** 2

        norm = math.sqrt(norm)

        if norm < 1e-12:
            continue

        if R[k][k] >= 0:
            alpha = -norm
        else:
            alpha = norm

        v = [0.0] * n

        for i in range(k, n):
            v[i] = R[i][k]

        v[k] -= alpha

        v_norm = 0.0

        for i in range(k, n):
            v_norm += v[i] ** 2

        if v_norm < 1e-24:
            continue

        for j in range(k, n):
            s = 0.0

            for i in range(k, n):
                s += v[i] * R[i][j]

            for i in range(k, n):
                R[i][j] -= 2.0 * v[i] * s / v_norm

        for i in range(n):
            s = 0.0

            for j in range(k, n):
                s += Q[i][j] * v[j]

            for j in range(k, n):
                Q[i][j] -= 2.0 * s * v[j] / v_norm

    return Q, R


def qr_algorithm(A, epsilon, max_iterations=10000):
    Ak = copy_matrix(A)

    for i in range(max_iterations):
        Q, R = qr_decomposition(Ak)
        Ak = multiply(R, Q)

        total = 0.0

        for i in range(len(Ak)):
            for j in range(i):
                total += Ak[i][j] ** 2

        if math.sqrt(total) < epsilon:
            break

    eigenvalues = []
    i = 0
    n = len(Ak)

    while i < n:
        if i == n - 1 or abs(Ak[i + 1][i]) < epsilon:
            eigenvalues.append(complex(Ak[i][i], 0))
            i += 1
        else:
            a = Ak[i][i]
            b = Ak[i][i + 1]
            c = Ak[i + 1][i]
            d = Ak[i + 1][i + 1]

            trace = a + d
            determinant = a * d - b * c
            discriminant = trace ** 2 - 4 * determinant

            if discriminant >= 0:
                root = math.sqrt(discriminant)
                eigenvalues.append(complex((trace + root) / 2, 0))
                eigenvalues.append(complex((trace - root) / 2, 0))
            else:
                root = math.sqrt(-discriminant)
                eigenvalues.append(complex(trace / 2, root / 2))
                eigenvalues.append(complex(trace / 2, -root / 2))

            i += 2

    return eigenvalues


def main():
    epsilon, A = read_input("5.txt")

    print("Точность:", epsilon)

    print()
    print("Исходная матрица A:")
    print_matrix(A)

    Q, R = qr_decomposition(A)

    print()
    print("Матрица Q:")
    print_matrix(Q)

    print()
    print("Матрица R:")
    print_matrix(R)

    QR = multiply(Q, R)

    print()
    print("Q * R:")
    print_matrix(QR)

    eigenvalues = qr_algorithm(A, epsilon)

    print()
    print("Собственные значения A:")

    for value in eigenvalues:
        print(value)

    _, B = read_input("complex.txt", False)

    print()
    print("Матрица B:")
    print_matrix(B)

    import numpy as np

    values = np.linalg.eigvals(B)

    print()
    print("Собственные значения B:")

    for value in values:
        print(value)

    if len(values) == 2:
        if (abs(values[0].real - values[1].real) < 1e-10
                and abs(values[0].imag + values[1].imag) < 1e-10
                and abs(values[0].imag) > 1e-10):
            print("Комплексно-сопряженная пара есть")
        else:
            print("Комплексно-сопряженная пара не найдена")


main()