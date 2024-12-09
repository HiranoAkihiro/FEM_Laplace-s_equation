nnode = 0
ndof = 0
node = []
with open('node.dat', 'r') as f:
    nnode, ndof =map(int, f.readline().split())

    for line in f:
        x, y = map(float, line.split())
        node.append([x,y])

id = []
id_e = 0
id_i = 0
for i in range(nnode):
    if node[i][1] == 2.0:
        if node[i][0] == 0.0:
            id_i = i+1
        elif node[i][0] == 1.0:
            id_e = i+1
        else :
            id.append(i+1)

force = 100/(2*(len(id)+2)-2)
with open('load.dat', 'w') as f:
    f.write(f"{len(id)+2} 3\n")
    f.write(f"{id_i}, 2, {force}\n")
    f.write(f"{id_e}, 2, {force}\n")
    for i in range(len(id)):
        f.write(f"{id[i]}, 2, {2.0*force}\n")

id_x = []
id_y = []
for i in range(nnode):
    if node[i][1] == 0.0:
        id_y.append(i+1)

for i in range(nnode):
    if node[i][0] == 0.0:
        id_x.append(i+1)

with open('bc.dat', 'w') as f:
    f.write(f"{len(id_x)+len(id_y)}, 3\n")
    for i in range(len(id_x)):
        f.write(f"{id_x[i]}, 1, 0.0\n")

    for i in range(len(id_y)):
        f.write(f"{id_y[i]}, 2, 0.0\n")

print(node)
print(id)
print(id_e, id_i)