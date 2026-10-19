import math


def read_matrix(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    A = []

    for line in lines:
        row = []

        for value in line.split():
            row.append(float(value))

        if len(row) > 0:
            A.append(row)

    return A


def copy_matrix(A):
    B = []

    for i in range(len(A)):
        row = []

        for j in range(len(A[i])):
            row.append(A[i][j])

        B.append(row)

    return B


def print_matrix(A):
    for row in A:
        for value in row:
            print(round(value, 6), end=" ")

        print()


def dot_product(a, b):
    s = 0.0

    for i in range(len(a)):
        s += a[i] * b[i]

    return s


def qr_decomposition(A):
    n = len(A)

    Q = []
    R = []

    for i in range(n):
        Q.append([0.0] * n)
        R.append([0.0] * n)

    for j in range(n):
        v = []

        for i in range(n):
            v.append(A[i][j])

        for k in range(j):
            q = []

            for i in range(n):
                q.append(Q[i][k])

            R[k][j] = dot_product(q, v)

            for i in range(n):
                v[i] -= R[k][j] * Q[i][k]

        norm = math.sqrt(dot_product(v, v))

        if norm < 1e-12:
            raise ValueError("Столбцы матрицы линейно зависимы")

        R[j][j] = norm

        for i in range(n):
            Q[i][j] = v[i] / norm

    return Q, R


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


def qr_algorithm(A, epsilon=1e-10, max_iterations=1000):
    Ak = copy_matrix(A)

    for iteration in range(max_iterations):
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
            discriminant = trace * trace - 4 * determinant

            if discriminant >= 0:
                root = math.sqrt(discriminant)
                eigenvalues.append(complex((trace + root) / 2, 0))
                eigenvalues.append(complex((trace - root) / 2, 0))
            else:
                root = math.sqrt(-discriminant)
                eigenvalues.append(
                    complex(trace / 2, root / 2)
                )
                eigenvalues.append(
                    complex(trace / 2, -root / 2)
                )

            i += 2

    return eigenvalues, Ak


def main():
    A = read_matrix("5.txt")

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
    print("Произведение Q * R:")
    print_matrix(QR)

    print()
    print("Собственные значения:")
    eigenvalues, Ak = qr_algorithm(A)

    for value in eigenvalues:
        print(value)


main()