import numpy as np

def nonzero(xs):
    return np.nonzero(xs)[0]

def nonfactor(xs):
    return nonzero(np.mod(xs[1:], xs[:-1]))

def is_one(x):
    return abs(x - 1) < 1e-09

def bezout(a, b):
    r0, r1 = a, b
    s0, s1 = 1, 0
    t0, t1 = 0, 1
    while r1 != 0:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1 
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if r0 < 0:
        r0, s0, t0 = -r0, -s0, -t0
    assert r0 == np.gcd(a, b), f"failed to compute gcd of {a} and {b}"
    assert r0 == a * s0 + b * t0, "failed to satisfy identity"
    return r0, s0, t0

# Apply row operations to get A[other_row, pivot_col] = 0
# Assumes A[pivot_row, pivot_col] does not divide A[other_row, pivot_col]
# Modifies: A, S
def eliminate_row_nonfactor(A, S, pivot_row, pivot_col, other_row):
    i, j, k = pivot_row, pivot_col, other_row
    beta, sig, tau = bezout(A[i,j], A[k,j])
    alpha, gamma = A[i,j] // beta, A[k,j] // beta
    L = np.array([[sig, tau], [-gamma, alpha]])
    A[[i, k]] = L @ A[[i, k]]
    S[[i, k]] = L @ S[[i, k]]
    return tau, gamma

# Apply row operations to get A[other_row, pivot_col] = 0
# Modifies: A, S
def eliminate_row(A, S, pivot_row, pivot_col, other_row):
    i, j, k = pivot_row, pivot_col, other_row
    if A[k,j] % A[i,j] == 0:
        q = A[k,j] // A[i,j]
        A[k] -= q * A[i]
        S[k] -= q * S[i]
    else:
        eliminate_row_nonfactor(A, S, i, j, k)

# Apply col operations to get A[pivot_row, other_col] = 0
# Modifies: A, T
def eliminate_col(A, T, pivot_row, pivot_col, other_col):
    eliminate_row(A.T, T.T, pivot_col, pivot_row, other_col)

# Wikipedia algorithm
# Could increase efficiency by not searching parts of the matrix we already know are zero
def smith_normal_form(matrix):
    A = np.copy(matrix)
    m, n = A.shape
    S, T = np.identity(m), np.identity(n)
    js = []
    j = -1
    for t in range(m):
        # Choose pivot (or break if there are none left)
        j += 1
        while j < n and len(ixs := nonzero(A[:,j])) == 0:
            j += 1
        if j == n: break
        if A[t,j] == 0:
            A[[t, ixs[0]]] = A[[ixs[0], t]]
            S[[t, ixs[0]]] = S[[ixs[0], t]]
        js.append(j)
        
        while True:
            # Eliminate entries in column
            for k in nonzero(A[:,j]):
                if k != t:
                    eliminate_row(A, S, t, j, k)
      
            # Eliminate entries in row
            for k in nonzero(A[t]):
                if k != j:
                    eliminate_col(A, T, t, j, k)
      
            # Repeat if necessary (I wish Python had goto)
            if len(nonzero(A[:,j])) == 1: break 

    r = len(js)

    # Move empty columns to the end (should be diagonal after this)
    for i in range(r):
        if not np.any(A[:,i]):
            A[:,[i,js[i]]] = A[:,[js[i],i]]
            T[:,[i,js[i]]] = T[:,[js[i],i]]

    # Ensure that entries are ordered by divisibility
    while len(ixs := nonfactor(np.diag(A)[:r])) > 0:
        i = ixs[0]
        A[:,i] += A[:,i+1]
        T[:,i] += T[:,i+1]
        tau, gamma = eliminate_row_nonfactor(A, S, i, i, i+1)
        A[:,i+1] -= tau * gamma * A[:,i]
        T[:,i+1] -= tau * gamma * T[:,i]

    # Make entries positive
    for i in range(r):
        if A[i, i] < 0:
            A[i] *= -1
            S[i] *= -1

    nonz = np.argwhere(A)
    correct_nonz = np.array([range(r), range(r)]).T
    detS, detT = np.linalg.det(S), np.linalg.det(T)
    assert np.all(A == S @ matrix @ T), "failed to satisfy identity"
    assert np.all(nonz == correct_nonz), "matrix A has nonzero entries in the wrong places"
    assert is_one(abs(detS)), f"matrix S failed to be invertible (determinant {detS})"
    assert is_one(abs(detT)), f"matrix T failed to be invertible (determinant {detT})"
    return np.diag(A)[:r]