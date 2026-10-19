def print_matrix(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            print(f"{A[i][j]:.6f}", end=" ")
        print()
    print()


def print_vector(x):
    for i in range(len(x)):
        print(f"x{i + 1} = {x[i]:.6f}")
    print()


def matrix_multiply(A, B):
    C = []

    for i in range(len(A)):
        row = []

        for j in range(len(B[0])):
            s = 0.0

            for k in range(len(B)):
                s += A[i][k] * B[k][j]

            row.append(s)

        C.append(row)

    return C


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


def inverse_from_lu(L, U, P):
    n = len(L)

    A_inv = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(0.0)

        A_inv.append(row)

    for j in range(n):
        e = []

        for i in range(n):
            if i == j:
                e.append(1.0)
            else:
                e.append(0.0)

        x = solve_lu(L, U, P, e)

        for i in range(n):
            A_inv[i][j] = x[i]

    return A_inv


def determinant_from_lu(U, swaps):
    det = 1.0

    for i in range(len(U)):
        det *= U[i][i]

    det *= (-1) ** swaps

    return det


def main():
    with open("1.txt", "r", encoding="utf-8") as file:
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

    L, U, P, swaps = lu_decomposition(A)

    LU = matrix_multiply(L, U)
    PA = matrix_multiply(P, A)
    x = solve_lu(L, U, P, b)
    A_inv = inverse_from_lu(L, U, P)
    det_A = determinant_from_lu(U, swaps)
    AA_inv = matrix_multiply(A, A_inv)

    print("L:")
    print_matrix(L)

    print("U:")
    print_matrix(U)

    print("L * U:")
    print_matrix(LU)

    print("P * A:")
    print_matrix(PA)

    print("x:")
    print_vector(x)

    print("A^(-1):")
    print_matrix(A_inv)

    print(f"det(A) = {det_A:.6f}")
    print()

    print("A * A^(-1):")
    print_matrix(AA_inv)


main()