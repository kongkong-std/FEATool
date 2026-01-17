import argparse
import numpy as np
import os
from pathlib import Path

def vertex_coordinate_process_dist(num_proc, file_path):
    vtxdist = [0]
    vertex_coor = []
    with open(file_path, 'r') as f:
        nv = int(f.readline().strip())

        q, r = divmod(nv, num_proc)
        local_size = [q] * num_proc
        for index in range(r):
            local_size[index] += 1

        for size in local_size:
            vtxdist.append(vtxdist[-1] + size)

        for _ in range(nv):
            line = f.readline()
            if not line:
                break
            parts = line.strip().split()
            idx = int(parts[0])
            vertex_type = int(parts[1])
            value = [float(x) for x in parts[2:]]
            vertex_coor.append([idx, vertex_type, value])

    # distributed vertex coordinate data
    vertex_coor_dist = []
    for rank in range(num_proc):
        index_start = vtxdist[rank]
        index_end = vtxdist[rank + 1]
        local_data = vertex_coor[index_start : index_end]
        vertex_coor_dist.append(local_data)

    return np.array(vtxdist, dtype=np.int32), np.array(vertex_coor_dist, dtype=object)

def write_vertex_coordinate_process_dist(file_path, vertex_coor_dist):
    num_proc = len(vertex_coor_dist)

    base_path = Path(file_path)
    local_path = base_path / f"np-{num_proc}"
    local_path.mkdir(parents=True, exist_ok=True)

    for index in range(num_proc):
        local_file = local_path / f"vertex_coor-{index}.txt"
        with open(local_file, 'w') as f:
            local_nv = len(vertex_coor_dist[index])
            f.write(f"{local_nv}\n")
            for val in vertex_coor_dist[index]:
                idx = int(val[0])
                vertex_type = int(val[1])
                x, y, z, nx, ny, nz = val[2]
                f.write(f"{idx}\t{vertex_type}\t")
                f.write(f"{x:021.16e}\t{y:021.16e}\t{z:021.16e}\t")
                f.write(f"{nx:021.16e}\t{ny:021.16e}\t{nz:021.16e}\n")

def adjacent_list_process_dist(vtxdist, file_path):
    num_proc = len(vtxdist) - 1
    vtx_ptr = []
    adj_vtx_idx = []

    with open(file_path, 'r') as f:
        nv = int(f.readline().strip())    # number of vertices

        for _ in range(nv + 1):
            line = f.readline()
            if not line:
                break
            vtx_ptr.append(int(line.strip()))

        nadj_v = vtx_ptr[nv] - vtx_ptr[0]    # number of adjacent vertices
        for _ in range(nadj_v):
            line = f.readline()
            if not line:
                break
            adj_vtx_idx.append(int(line.strip()))

    vtx_ptr_dist = []    # distributed vertex pointer
    adj_vtx_idx_dist = []    # distributed adjacent vertex list
    for rank in range(num_proc):
        index_start = vtxdist[rank]
        index_end = vtxdist[rank + 1]
        vtx_ptr_start = vtx_ptr[index_start]
        vtx_ptr_end = vtx_ptr[index_end]
        
        local_adj_vtx_idx = adj_vtx_idx[vtx_ptr_start : vtx_ptr_end]
        adj_vtx_idx_dist.append(local_adj_vtx_idx)

        local_vtx_ptr = [ptr - vtx_ptr_start for ptr in vtx_ptr[index_start : index_end + 1]]
        vtx_ptr_dist.append(local_vtx_ptr)

    return np.array(vtx_ptr_dist, dtype=object), np.array(adj_vtx_idx_dist, dtype=object)

def write_adjacent_list_process_dist(file_path, vtxdist, vtx_ptr_dist, adj_vtx_idx_dist):
    num_proc = len(vtx_ptr_dist)

    base_path = Path(file_path)
    local_path = base_path / f"np-{num_proc}"
    local_path.mkdir(parents=True, exist_ok=True)

    for index in range(num_proc):
        local_file = local_path / f"adjacent_list-{index}.txt"
        with open(local_file, 'w') as f:
            local_nv = len(vtx_ptr_dist[index]) - 1
            f.write(f"{local_nv}\n")
            for val in list(range(vtxdist[index], vtxdist[index + 1])):
                f.write(f"{val}\n")
            for val in vtx_ptr_dist[index]:
                f.write(f"{val}\n") 
            for val in adj_vtx_idx_dist[index]:
                f.write(f"{val}\n")

def matrix_process_dist(vertex_coor_dist, file_path):
    num_proc = len(vertex_coor_dist)

    # ---- Read CSR matrix ----
    row_ptr = []
    col_idx = []
    val = []

    with open(file_path, 'r') as f:
        parts = f.readline().strip().split()
        m, n, nnz = map(int, parts)

        for _ in range(m + 1):
            row_ptr.append(int(f.readline().strip()))

        for _ in range(nnz):
            col_idx.append(int(f.readline().strip()))

        for _ in range(nnz):
            val.append(float(f.readline().strip()))

    # ---- Build dof distribution table ----
    dofdist = [0]
    for rank in range(num_proc):
        local_dof = 0
        for v in vertex_coor_dist[rank]:
            local_dof += 6 if v[1] == 0 else 3
        dofdist.append(dofdist[-1] + local_dof)

    # ---- Distributed CSR arrays ----
    row_idx_dist = []
    row_ptr_dist = []
    col_idx_dist = []
    val_dist = []

    # ---- Process each rank ----
    for rank in range(num_proc):
        global_start_dof = dofdist[rank]
        global_end_dof   = dofdist[rank + 1] - 1
        local_dof = dofdist[rank + 1] - dofdist[rank]

        # row_ptr slice indices
        rp_row_start = global_start_dof
        rp_row_end   = global_end_dof + 1  # inclusive row, so +1

        rp_start = row_ptr[rp_row_start]
        rp_end   = row_ptr[rp_row_end]

        local_row_idx = list(range(global_start_dof, global_end_dof + 1))
        row_idx_dist.append(local_row_idx)

        # build local row_ptr (size = local_dof + 1)
        local_row_ptr = [
            row_ptr[rp_row_start + i] - rp_start
            for i in range(local_dof + 1)
        ]
        row_ptr_dist.append(local_row_ptr)

        # col_idx and val segments
        local_col_idx = col_idx[rp_start : rp_end]
        col_idx_dist.append(local_col_idx)

        local_val = val[rp_start : rp_end]
        val_dist.append(local_val)

    return (
        np.array(row_idx_dist, dtype=object),
        np.array(row_ptr_dist, dtype=object),
        np.array(col_idx_dist, dtype=object),
        np.array(val_dist, dtype=object)
    )

def write_matrix_process_dist(file_path, row_idx_dist, row_ptr_dist, col_idx_dist, val_dist):
    num_proc = len(row_ptr_dist)

    base_path = Path(file_path)
    local_path = base_path / f"np-{num_proc}"
    local_path.mkdir(parents=True, exist_ok=True)

    for index in range(num_proc):
        local_file = local_path / f"A-{index}.txt"
        local_row_idx = row_idx_dist[index]
        local_row_ptr = row_ptr_dist[index]
        local_col_idx = col_idx_dist[index]
        local_val = val_dist[index]
        with open(local_file, 'w') as f:
            local_m = len(local_row_ptr) - 1
            local_nnz = local_row_ptr[local_m] - local_row_ptr[0]
            f.write(f"{local_m}\t{local_m}\t{local_nnz}\n")
            for tmp_val in local_row_idx:
                f.write(f"{tmp_val}\n")
            for tmp_val in local_row_ptr:
                f.write(f"{tmp_val}\n")
            for tmp_val in local_col_idx:
                f.write(f"{tmp_val}\n")
            for tmp_val in local_val:
                f.write(f"{tmp_val:021.16e}\n")

def vector_process_dist(vertex_coor_dist, file_path):
    num_proc = len(vertex_coor_dist)

    # --------------------------
    # read full vector
    # --------------------------
    val = []
    with open(file_path, 'r') as f:
        m = int(f.readline().strip())
        for _ in range(m):
            val.append(float(f.readline().strip()))

    # --------------------------
    # build dofdist (similar to vtxdist)
    # --------------------------
    dofdist = [0]
    for rank in range(num_proc):
        local_dof = 0
        for v in vertex_coor_dist[rank]:
            if v[1] == 0:
                local_dof += 6
            elif v[1] == 1:
                local_dof += 3
            else:
                raise ValueError("Unknown vertex type")

        dofdist.append(dofdist[-1] + local_dof)

    # --------------------------
    # distribute vector
    # --------------------------
    idx_dist = []
    val_dist = []

    for rank in range(num_proc):
        start = dofdist[rank]
        end   = dofdist[rank + 1]  # exclusive

        local_idx = list(range(start, end))
        idx_dist.append(local_idx)

        local_vec = val[start : end]
        val_dist.append(local_vec)

    return np.array(idx_dist, dtype=object), np.array(val_dist, dtype=object)

def write_vector_process_dist(file_path, idx_dist, vec_dist):
    num_proc = len(vec_dist)

    base_path = Path(file_path)
    local_path = base_path / f"np-{num_proc}"
    local_path.mkdir(parents=True, exist_ok=True)

    for index in range(num_proc):
        local_file = local_path / f"b-{index}.txt"
        local_idx = idx_dist[index]
        local_val = vec_dist[index]
        with open(local_file, 'w') as f:
            local_m = len(local_val)
            f.write(f"{local_m}\n")
            for tmp_idx, tmp_val in zip(local_idx, local_val):
                f.write(f"{tmp_idx}\t{tmp_val:021.16e}\n")
            #for tmp_val in local_idx:
            #    f.write(f"{tmp_val}\n")
            #for tmp_val in local_val:
            #    f.write(f"{tmp_val:021.16e}\n")

def parse_args():
    parser = argparse.ArgumentParser(
        description="parse mesh and linear system data"
    )

    parser.add_argument(
        "-n", "--num_proc",
        type=int,
        required=True,
        help="Number of processors"
    )

    parser.add_argument(
        "-v", "--vertex_coor",
        type=str,
        required=True,
        help="Path to mesh vertex coordinate file"
    )

    parser.add_argument(
        "--output_vertex_coor",
        type=str,
        required=True,
        help="Output path for vertex coordinate"
    )

    parser.add_argument(
        "-a", "--adjacent_list",
        type=str,
        required=True,
        help="Path to mesh adjacent list file"
    )

    parser.add_argument(
        "--output_adjacent_list",
        type=str,
        required=True,
        help="Output path for adjacent list"
    )

    parser.add_argument(
        "--mat_file",
        type=str,
        required=True,
        help="Path to matrix file"
    )

    parser.add_argument(
        "--output_mat_file",
        type=str,
        required=True,
        help="Output path for matrix file"
    )

    parser.add_argument(
        "--vec_file",
        type=str,
        required=True,
        help="Path to vector file"
    )

    parser.add_argument(
        "--output_vec_file",
        type=str,
        required=True,
        help="Output path for vector file"
    )

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    print("Number of processors:", args.num_proc)
    print("Mesh coordinate file:", args.vertex_coor)
    print("Output path for vertex coordinate:", args.output_vertex_coor)
    print("Mesh adjacent list file:", args.adjacent_list)
    print("Output path for adjacent list:", args.output_adjacent_list)
    print("Matrix file:", args.mat_file)
    print("Output path for matrix file:", args.output_mat_file)
    print("Vector file:", args.vec_file)
    print("Output path for vector file:", args.output_vec_file)

    # distributed vertex array and vertex coordinate data
    vtxdist, vertex_coor_dist = vertex_coordinate_process_dist(args.num_proc, args.vertex_coor)

    # output distributed data to file
    write_vertex_coordinate_process_dist(args.output_vertex_coor, vertex_coor_dist)

    # distributed mesh adjacent list data
    vtx_ptr_dist, adj_vtx_idx_dist = adjacent_list_process_dist(vtxdist, args.adjacent_list)

    # output distributed adjacent list data
    write_adjacent_list_process_dist(args.output_adjacent_list, vtxdist, vtx_ptr_dist, adj_vtx_idx_dist)

    # distributed matrix data
    mat_row_idx_dist, mat_row_ptr_dist, mat_col_idx_dist, mat_val_dist = matrix_process_dist(vertex_coor_dist, args.mat_file)

    # output distributed matrix data
    write_matrix_process_dist(args.output_mat_file, mat_row_idx_dist, mat_row_ptr_dist, mat_col_idx_dist, mat_val_dist)

    # distributed vector data
    vec_idx_dist, vec_val_dist = vector_process_dist(vertex_coor_dist, args.vec_file)

    # ouput distributed vector data
    write_vector_process_dist(args.output_vec_file, vec_idx_dist, vec_val_dist)

# usage
# python parfile-process.py --num_proc <number of processors> \
#     --vertex_coor </path/to/vertex/coordinate/file> \
#     --output_vertex_coor </output/path/to/vertex/coordinate/file> \
#     --adjacent_list </path/to/adjacent/list/file> \
#     --output_adjacent_list </output/path/to/adjacent/list/file> \
#     --mat_file </path/to/matrix/file> \
#     --output_mat_file </output/path/to/matrix/file> \
#     --vec_file </path/to/vector/file> \
#     --output_vec_file </output/path/to/vector/file>
#
