import numpy as np
from smith import smith_normal_form
from spaces import sphere, torus, klein_bottle

matrix = np.array([[2, 4, 4], [-6, 6, 12], [10, 4, 16]])
v = smith_normal_form(matrix)
assert np.all(v == np.array([2, 2, 156])), "smith_normal_form fails test 1"

matrix = np.array([[-6, 111, -36, 6], [5, -672, 210, 74], [0, -255, 81, 24], [-7, 255, -81, -10]])
v = smith_normal_form(matrix)
assert np.all(v == np.array([1, 3, 21])), "smith_normal_form fails test 2"

matrix = np.array([
    3 * [-1] + 6 * [0], 
    [1, 0, 0, -1, -1, 0, 0, 0, 0], 
    [0, 1, 0, 1, 0, -1, -1, 0, 0],
    [0, 0, 1, 0, 1, 1, 0, -1, 0],
    6 * [0] + [1, 1, 0], 
    8 * [0] + [-1],
    8 * [0] + [1]
])
v = smith_normal_form(matrix)
assert np.all(v == 1), "smith_normal_form fails test 3"

H = [ G.pretty(False) for G in sphere(3).homology() ]
assert H == ["Z", "0", "0", "Z"], "homology fails for 3-sphere"

H = [ G.pretty(False) for G in torus().homology() ]
assert H == ["Z", "Z^2", "Z"], "homology fails for Torus"

H = [ G.pretty(False) for G in klein_bottle().homology() ]
assert H == ["Z", "Z + Z/2Z", "0"], "homology fails for Klein bottle"

print("All tests passed!")
