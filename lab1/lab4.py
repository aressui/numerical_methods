import math

def read_input(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    epsilon = float(lines[0])

    A = []

    for line in lines[1:]:
        row = []

        for value in line.split():
            row.append(float(value))

        A.append(row)

    return epsilon, A


def max_offdiagonal(A):
    n = len(A)

    p = 0
    q = 1

    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j]) > abs(A[p][q]):
                p = i
                q = j

    return p, q


def offdiagonal_measure(A):
    n = len(A)

    s = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            s += A[i][j] ** 2

    return math.sqrt(s)


def jacobi(A, epsilon):
    n = len(A)

    A_k = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(A[i][j])

        A_k.append(row)

    V = []

    for i in range(n):
        row = []

        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                row.append(0.0)

        V.append(row)

    rotations = 0

    t = offdiagonal_measure(A_k)

    while t >= epsilon:
        p, q = max_offdiagonal(A_k)

        if A_k[p][p] == A_k[q][q]:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(
                2.0 * A_k[p][q] /
                (A_k[p][p] - A_k[q][q])
            )

        c = math.cos(phi)
        s = math.sin(phi)

        old_A = []

        for i in range(n):
            row = []

            for j in range(n):
                row.append(A_k[i][j])

            old_A.append(row)

        for i in range(n):
            if i != p and i != q:
                A_k[i][p] = (
                    c * old_A[i][p] +
                    s * old_A[i][q]
                )

                A_k[p][i] = A_k[i][p]

                A_k[i][q] = (
                    -s * old_A[i][p] +
                    c * old_A[i][q]
                )

                A_k[q][i] = A_k[i][q]

        A_k[p][p] = (
            c * c * old_A[p][p] +
            2.0 * c * s * old_A[p][q] +
            s * s * old_A[q][q]
        )

        A_k[q][q] = (
            s * s * old_A[p][p] -
            2.0 * c * s * old_A[p][q] +
            c * c * old_A[q][q]
        )

        A_k[p][q] = 0.0
        A_k[q][p] = 0.0

        old_V = []

        for i in range(n):
            row = []

            for j in range(n):
                row.append(V[i][j])

            old_V.append(row)

        for i in range(n):
            V[i][p] = (
                c * old_V[i][p] +
                s * old_V[i][q]
            )

            V[i][q] = (
                -s * old_V[i][p] +
                c * old_V[i][q]
            )

        rotations += 1

        t = offdiagonal_measure(A_k)

    eigenvalues = []

    for i in range(n):
        eigenvalues.append(A_k[i][i])

    return eigenvalues, V, rotations


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


def print_matrix(A):
    for row in A:
        print(*row)


def main():
    epsilon, A = read_input("4.txt")

    eigenvalues, V, rotations = jacobi(A, epsilon)

    print("Заданная точность:", epsilon)

    print()
    print("Собственные значения:")

    for i in range(len(eigenvalues)):
        print("lambda_", i + 1, "=", eigenvalues[i])

    print()
    print("Собственные векторы:")

    for j in range(len(V)):
        print("v_", j + 1, "=", end=" ")

        for i in range(len(V)):
            print(V[i][j], end=" ")

        print()

    print()
    print("Матрица собственных векторов V:")
    print_matrix(V)

    Lambda = []

    for i in range(len(A)):
        row = []

        for j in range(len(A)):
            if i == j:
                row.append(eigenvalues[i])
            else:
                row.append(0.0)

        Lambda.append(row)

    AV = multiply(A, V)
    VL = multiply(V, Lambda)

    print()
    print("A * V:")
    print_matrix(AV)

    print()
    print("V * Lambda:")
    print_matrix(VL)

    difference = 0.0

    for i in range(len(A)):
        for j in range(len(A)):
            value = abs(AV[i][j] - VL[i][j])

            if value > difference:
                difference = value

    print()
    print("Максимальная разность:", difference)

    print()
    print("Количество вращений:", rotations)


main()