def print_vector(x, name):
    print(name)

    for i in range(len(x)):
        print(f"{name}{i + 1} = {x[i]:.6f}")

    print()


def matrix_norm(matrix):
    q = 0.0

    for i in range(len(matrix)):
        s = 0.0

        for j in range(len(matrix)):
            s += abs(matrix[i][j])

        if s > q:
            q = s

    return q


def make_iteration_matrix(A, b):
    n = len(A)

    alpha = []
    beta = []

    for i in range(n):
        row = []

        for j in range(n):
            if i == j:
                row.append(0.0)
            else:
                row.append(-A[i][j] / A[i][i])

        alpha.append(row)
        beta.append(b[i] / A[i][i])

    return alpha, beta


def simple_iterations(alpha, beta, epsilon):
    n = len(beta)

    q = matrix_norm(alpha)

    x = beta[:]

    for k in range(10000):
        x_new = [0.0] * n

        for i in range(n):
            s = beta[i]

            for j in range(n):
                s += alpha[i][j] * x[j]

            x_new[i] = s

        difference = 0.0

        for i in range(n):
            value = abs(x_new[i] - x[i])

            if value > difference:
                difference = value

        error = q / (1 - q) * difference

        if error <= epsilon:
            return x_new, k + 1, q

        x = x_new

    return x, 10000, q


def seidel(alpha, beta, epsilon):
    n = len(beta)

    x = beta[:]

    G = [[0.0] * n for i in range(n)]

    for j in range(n):
        old = [0.0] * n
        old[j] = 1.0

        new = old[:]

        for i in range(n):
            s = 0.0

            for m in range(i):
                s += alpha[i][m] * new[m]

            for m in range(i + 1, n):
                s += alpha[i][m] * old[m]

            new[i] = s

        for i in range(n):
            G[i][j] = new[i]

    q = matrix_norm(G)

    for k in range(10000):
        x_old = x[:]

        for i in range(n):
            s = beta[i]

            for j in range(i):
                s += alpha[i][j] * x[j]

            for j in range(i + 1, n):
                s += alpha[i][j] * x_old[j]

            x[i] = s

        difference = 0.0

        for i in range(n):
            value = abs(x[i] - x_old[i])

            if value > difference:
                difference = value

        error = q / (1 - q) * difference

        if error <= epsilon:
            return x, k + 1, q

    return x, 10000, q


def lu_decomposition(A):
    n = len(A)

    U = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(A[i][j])

        U.append(row)

    L = []

    for i in range(n):
        row = []

        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                row.append(0.0)

        L.append(row)

    P = []

    for i in range(n):
        row = []

        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                row.append(0.0)

        P.append(row)

    swaps = 0

    for k in range(n - 1):
        p = k

        for i in range(k + 1, n):
            if abs(U[i][k]) > abs(U[p][k]):
                p = i

        if p != k:
            temp = U[k]
            U[k] = U[p]
            U[p] = temp

            temp = P[k]
            P[k] = P[p]
            P[p] = temp

            for j in range(k):
                temp = L[k][j]
                L[k][j] = L[p][j]
                L[p][j] = temp

            swaps += 1

        for i in range(k + 1, n):
            L[i][k] = U[i][k] / U[k][k]

            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

            U[i][k] = 0.0

    return L, U, P, swaps


def solve_lu(L, U, P, b):
    n = len(L)

    Pb = []

    for i in range(n):
        s = 0.0

        for j in range(n):
            s += P[i][j] * b[j]

        Pb.append(s)

    z = []

    for i in range(n):
        s = 0.0

        for j in range(i):
            s += L[i][j] * z[j]

        z.append(Pb[i] - s)

    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        s = 0.0

        for j in range(i + 1, n):
            s += U[i][j] * x[j]

        x[i] = (z[i] - s) / U[i][i]

    return x


def main():
    with open("3.txt", "r", encoding="utf-8") as file:
        data = file.read().split()

    p = 0

    n = int(data[p])
    p += 1

    A = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(float(data[p]))
            p += 1

        A.append(row)

    b = []

    for i in range(n):
        b.append(float(data[p]))
        p += 1

    epsilon = float(data[p])

    alpha, beta = make_iteration_matrix(A, b)

    L, U, P, swaps = lu_decomposition(A)
    x_exact = solve_lu(L, U, P, b)

    x_simple, k_simple, q_simple = simple_iterations(
        alpha, beta, epsilon
    )

    x_seidel, k_seidel, q_seidel = seidel(
        alpha, beta, epsilon
    )

    print(f"Заданная точность epsilon = {epsilon}")
    print()

    print("Метод простых итераций")
    print(f"Норма q = {q_simple:.6f}")
    print(f"Количество итераций = {k_simple}")
    print_vector(x_simple, "x")

    print("Метод Зейделя")
    print(f"Норма q = {q_seidel:.6f}")
    print(f"Количество итераций = {k_seidel}")
    print_vector(x_seidel, "x")

    print("Эталонное решение")
    print_vector(x_exact, "x")

    print("Сравнение методов")
    print(f"Метод простых итераций: {k_simple} итераций")
    print(f"Метод Зейделя: {k_seidel} итераций")


main()