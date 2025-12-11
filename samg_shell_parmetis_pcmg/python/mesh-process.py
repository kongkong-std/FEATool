import argparse
import numpy as np

def read_entity(file_path):
    # read entity index data
    entity = []
    ns = 0    # number of shell entities
    flag_shell = False
    cnt_shell = 0

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if "shell" in line:
                flag_shell = True
                ns = int(line.split()[1])
                continue

            if flag_shell and cnt_shell < ns:
                parts = line.split()
                entity.append(int(parts[0]))
                cnt_shell += 1
                if cnt_shell == ns:
                    break
            continue

    return np.array(entity, dtype=np.int32)

def write_adjacent_list(ele_triangle, ele_tetrahedron, file_path):
    # write element data to file in CSR format

    # number of vertices in mesh
    if ele_triangle.size > 0:
        max_tri = ele_triangle[:, 1:].max()
    else:
        max_tri = -1

    if ele_tetrahedron.size > 0:
        max_tet = ele_tetrahedron[:, 1:].max()
    else:
        max_tet = -1

    nv = int(max(max_tri, max_tet)) + 1
    neighbours = [set() for _ in range(nv)]

    # triangle element
    for elem in ele_triangle:
        nodes = elem[1:]
        for index_i in range(3):
            for index_j in range(3):
                if index_i != index_j:
                    neighbours[nodes[index_i]].add(nodes[index_j])

    # tetrahedron element
    for elem in ele_tetrahedron:
        nodes = elem[1:]
        for index_i in range(4):
            for index_j in range(4):
                if index_i != index_j:
                    neighbours[nodes[index_i]].add(nodes[index_j])

    xadj = np.zeros(nv + 1, dtype=np.int32)
    adjncy_list = []
    for index_i in range(nv):
        nbrs = sorted(neighbours[index_i])
        adjncy_list.extend(nbrs)
        xadj[index_i + 1] = xadj[index_i] + len(nbrs)

    adjncy_list = np.array(adjncy_list, dtype=np.int32)

    with open(file_path, 'w') as f:
        f.write(f"{nv}\n")
        for val in xadj:
            f.write(f"{val}\n")
        for val in adjncy_list:
            f.write(f"{val}\n")
    

def read_element(file_path, element_name, node_count):
    ele = []
    ele_entity = []
    ne = None
    ne_entity = None

    flag_section = False
    flag_ele = False
    flag_entity = False

    cnt_ele = 0
    cnt_entity = 0
            
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # 开始对应单元 section
            if f"{element_name} # type name" in line:
                flag_section = True
                continue

            if not flag_section:
                continue

            if "# number of elements" in line:
                ne = int(line.split()[0])
                continue

            if "# Elements" in line:
                flag_ele = True
                continue

            if flag_ele and cnt_ele < ne:
                parts = line.split()
                ele.append([int(p) for p in parts[:node_count]])
                cnt_ele += 1
                if cnt_ele == ne:
                    flag_ele = False
                continue

            if "# number of geometric entity indices" in line:
                ne_entity = int(line.split()[0])
                continue

            if "# Geometric entity indices" in line:
                flag_entity = True
                continue

            if flag_entity and cnt_entity < ne_entity:
                ele_entity.append(int(line.split()[0]))
                cnt_entity += 1
                if cnt_entity == ne_entity:
                    flag_entity = False
                continue

    # ---------------------------------------------------------
    #                ✓ 关键增强：处理不存在的单元类型
    # ---------------------------------------------------------
    if ne is None and ne_entity is None:
        # 单元不存在，返回空数组（shape = 0 × (node_count+1)）
        return np.zeros((0, node_count + 1), dtype=np.int32)

    # sanity check
    if ne != ne_entity:
        raise ValueError("Element count mismatch")

    ele = np.array(ele, dtype=np.int32)
    ele_entity = np.array(ele_entity, dtype=np.int32)

    return np.column_stack((ele_entity, ele))

def read_triangle_element(file_path):
    return read_element(file_path, "tri", 3)

def read_tetrahedron_element(file_path):
    return read_element(file_path, "tet", 4)

def write_vertex_coordinate(vertices, file_path):
    # write vertex coordinate to file path
    nv = vertices.shape[0]

    with open(file_path, 'w') as f:
        f.write(f"{nv}\n")
        for index, vertex in enumerate(vertices):
            vertextype = int(vertex[0])
            x, y, z = vertex[1], vertex[2], vertex[3]
            nx, ny, nz = vertex[4], vertex[5], vertex[6]
            f.write(f"{index}\t{vertextype}\t")
            f.write(f"{x:021.16e}\t{y:021.16e}\t{z:021.16e}\t")
            f.write(f"{nx:021.16e}\t{ny:021.16e}\t{nz:021.16e}\n")

def read_vertex_coordinate(file_path):
    # vertices coordinate data
    vertices = []
    nv = 0
    flag_vertex_coor = False
    cnt_vertex = 0

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue    # skip empty line

            if "# number of mesh vertices" in line:
                nv = int(line.split()[0])    # number of vertices

            if "# Mesh vertex coordinates" in line:
                flag_vertex_coor = True
                continue

            if flag_vertex_coor:
                parts = line.split()
                vertices.append([float(p) for p in parts])
                cnt_vertex += 1
                if cnt_vertex >= nv:
                    break

    return np.array(vertices, dtype=np.float64)

def filter_ele_triangle_entity(ele_triangle, entity):
    if ele_triangle.size == 0:
        return ele_triangle

    entity_set = set(entity.tolist())  # 转成 set 加速查找
    mask = np.array([row[0] in entity_set for row in ele_triangle], dtype=bool)
    return ele_triangle[mask]

def update_vertex_coordinate_type(vertex_coor, ele_triangle, ele_tetrahedron):

    nv = vertex_coor.shape[0]

    # 顶点类型初始化为 -1（没有出现）
    vertex_type = -np.ones(nv, dtype=np.int32)

    # --- 处理三角形：shell 顶点 type = 0 ---
    if ele_triangle.size > 0:
        tri_nodes = np.unique(ele_triangle[:, 1:].flatten())
        vertex_type[tri_nodes] = 0

    # --- 处理四面体：solid 顶点 type = 1 ---
    if ele_tetrahedron.size > 0:
        tet_nodes = np.unique(ele_tetrahedron[:, 1:].flatten())
        vertex_type[tet_nodes] = 1

    # 组合输出
    vertex_new = np.hstack((
        vertex_type.reshape(-1, 1),
        vertex_coor
    ))

    return vertex_new

def update_vertex_coordinate_normal(vertex_coor, ele_triangle):

    """
    vertex_coor: shape = (N, 4)
        每行结构 [type, x, y, z]
        type: 0 shell, 1 solid
    ele_triangle: shape = (M, 4)
        三角形单元顶点索引（0-based）
    """

    n_vertices = vertex_coor.shape[0]

    # 抽取 xyz 坐标
    coords = vertex_coor[:, 1:4]

    # 顶点法向累加器
    normal_sum = np.zeros((n_vertices, 3), dtype=float)

    # ==== 1. 面积加权计算 shell 顶点法向 ====
    for elem in ele_triangle:
        i1, i2, i3 = int(elem[1]), int(elem[2]), int(elem[3])

        r1 = coords[i1]
        r2 = coords[i2]
        r3 = coords[i3]

        # 三角形法向
        cross_vec = np.cross(r2 - r1, r3 - r1)
        area2 = np.linalg.norm(cross_vec)  # 实际面积 = area2/2

        if area2 == 0:
            continue

        n_elem = cross_vec / area2      # 单位法向（用未除2面积权重也 OK）
        weight = area2 * 0.5            # 面积 = |cross| / 2

        # 累加到三个顶点
        normal_sum[i1] += weight * n_elem
        normal_sum[i2] += weight * n_elem
        normal_sum[i3] += weight * n_elem

    # ==== 2. 归一化 + solid 顶点设为 0 ====
    normals = np.zeros_like(normal_sum)

    for i in range(n_vertices):
        if int(vertex_coor[i, 0]) == 0:  # shell 顶点
            norm = np.linalg.norm(normal_sum[i])
            if norm > 1e-14:
                normals[i] = normal_sum[i] / norm
            else:
                normals[i] = np.array([0.0, 0.0, 1.0])  # 极少退化情况
        else:
            normals[i] = np.array([0.0, 0.0, 0.0])  # solid 顶点 normal = 0

    # ==== 3. 合并，形成 [type x y z nx ny nz] ====
    vertex_updated = np.hstack((vertex_coor, normals))

    return vertex_updated

def parse_args():
    parser = argparse.ArgumentParser(
        description="parse mesh file and export data"
    )

    parser.add_argument(
        "-m", "--mesh",
        type=str,
        required=True,
        help="Path to mesh file"
    )

    parser.add_argument(
        "-v", "--vertex_coor",
        type=str,
        required=True,
        help="Output path for vertices coordinates"
    )

    parser.add_argument(
        "-a", "--adjacent_list",
        type=str,
        required=True,
        help="Output path for adjacent list CSR format"
    )

    parser.add_argument(
        "-e", "--entity",
        type=str,
        required=True,
        help="Path to entity file contains shell and solid entity"
    )

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    print("Mesh file:", args.mesh)
    print("Vertex coordinate output:", args.vertex_coor)
    print("Adjacent list output:", args.adjacent_list)
    print("Entity file:", args.entity)

    vertices_coor = read_vertex_coordinate(args.mesh)
    ele_triangle = read_triangle_element(args.mesh)
    ele_tetrahedron = read_tetrahedron_element(args.mesh)
    entities = read_entity(args.entity)

    ele_triangle = filter_ele_triangle_entity(ele_triangle, entities)    # filter triangle element with shell entities
    
    # updating vertices coordinate with type, type 0 in shell, type 1 in solid
    # computing vertices normal vector with element weighted average
    vertices_coor = update_vertex_coordinate_type(vertices_coor, ele_triangle, ele_tetrahedron)
    vertices_coor = update_vertex_coordinate_normal(vertices_coor, ele_triangle)

    write_vertex_coordinate(vertices_coor, args.vertex_coor)
    write_adjacent_list(ele_triangle, ele_tetrahedron, args.adjacent_list)
    print("Number of vertices:",len(vertices_coor))
    print("Number of triangle elements:", ele_triangle.shape[0])
    print("Number of tetrahedron elements:", ele_tetrahedron.shape[0])

# usage
# python mesh-process.py --mesh </path/to/mesh/file> \
# --vertex_coor </output/path/to/vertex/coordinate> \
# --adjacent_list </output/path/to/adjacent/list/csr/format>
# --entity </path/to/entity/file>
#
