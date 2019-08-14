# Computes the forward number-theoretic transform of the given vector,
# with respect to the given primitive nth root of unity under the given modulus.
# The length of the vector must be a power of 2.

def transform_fast(vector, root, mod):
    transform_vec = vector
    n = len(transform_vec)
    levels = n.bit_length() - 1
    if 1 << levels != n:
        raise ValueError("Length is not a power of 2")

    powtable = []
    temp = 1
    for i in range(n // 2):
        powtable.append(temp)
        temp = temp * root % mod

    # print("powtable[4]=",str(powtable[4]))
    def reverse(x, bits):
        y = 0
        for i in range(bits):
            y = (y << 1) | (x & 1)
            x >>= 1
        return y

    for i in range(n):
        j = reverse(i, levels)
        if j > i:
            transform_vec[i], transform_vec[j] = transform_vec[j], transform_vec[i]
    size = 2
    while size <= n:
        halfsize = size // 2
        tablestep = n // size
        for i in range(0, n, size):
            k = 0
            for j in range(i, i + halfsize):
                l = j + halfsize
                left = transform_vec[j]
                right = transform_vec[l] * powtable[k]
                transform_vec[j] = (left + right) % mod
                transform_vec[l] = (left - right) % mod
                k += tablestep
        size *= 2
    return (transform_vec)


# Compute the inverse of ntt

def inverse_transform(vector, inv_root, inv_L, mod):
    outvec = transform_fast(vector, inv_root, mod)
    return [(val * inv_L % mod) for val in outvec]


# Compute the convolution using ntt

def convolve_ntt(vec1, vec2, root, inv_root, inv_L, mod):
    temp1 = transform_fast(vec1, root, mod)
    temp2 = transform_fast(vec2, root, mod)
    temp3 = [(x * y % mod) for (x, y) in zip(temp1, temp2)]
    conv = inverse_transform(temp3, inv_root, inv_L, mod)
    return conv