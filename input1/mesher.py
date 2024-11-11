import sys
import numpy as np

# Constants
kint = 4
kdouble = 8

# Command line arguments
if len(sys.argv) == 3:
    EX = int(sys.argv[1])
    EY = int(sys.argv[2])
else:
    EX, EY = 1, 1

LX = 1.0
LY = 1.0
# Calculate dimensions
NELEM = EX * EY
NNODE = (EX+1) * (EY+1)
NX, NY = EX + 1, EY + 1

print(f"NNODE   : {NNODE:8d}")
print(f"NELEM   : {NELEM:8d}")
print(f"Length X: {LX:12.5e}")
print(f"Length Y: {LY:12.5e}")
print(f"Num of elem X: {EX:12d}")
print(f"Num of elem Y: {EY:12d}")

TX, TY = LX/EX, LY/EY

elemNode = np.zeros((3, NELEM*2), dtype=int)
coord = np.zeros((2, NNODE))

# Generate coordinates
in_ = 0
for j in range(NY):
    for i in range(NX):
        coord[:, in_] = [i*TX, j*TY]
        in_ += 1

# Generate element nodes
in_ = 0
for j in range(1,EY+1,1):
    for i in range(1,EX+1,1):
        elemNode[:, in_] = [
            (j-1)*NX + i,
            j*NX + i,
            j*NX + i + 1
        ]
        elemNode[:, in_+1] = [
            (j-1)*NX + i,
            (j-1)*NX + i + 1,
            j*NX + i + 1
        ]
        in_ += 2

# Write node data
with open("node.dat", "w") as f:
    f.write(f"{NNODE} 2\n")
    for i in range(NNODE):
        f.write(f"{coord[0,i]:12.5e},{coord[1,i]:12.5e}\n")

# Write element data
with open("elem.dat", "w") as f:
    f.write(f"{NELEM*2} 3\n")
    for i in range(NELEM*2):
        f.write(",".join(map(str, elemNode[:, i])) + "\n")

# Write boundary condition data
with open("bc.dat", "w") as f:
    f.write(f"{2} 3\n")
    f.write(f"{NX}, 1, 0.0\n")
    f.write(f"{NX*(NX-1)+1}, 1, 1.0\n")

# with open("bc.dat", "w") as f:
#     ip1 = sum(1 for k in range(NZ) for j in range(NY) for i in range(NX) if i == 0)
#     f.write(f"{3*ip1} 3\n")
#     in_ = 1
#     for j in range(NY):
#         for i in range(NX):
#             if i == 0:
#                 f.write(f"{in_}, 1, 0.0\n")
#                 f.write(f"{in_}, 2, 0.0\n")
#                 f.write(f"{in_}, 3, 0.0\n")
#             in_ += 1

# Write load data