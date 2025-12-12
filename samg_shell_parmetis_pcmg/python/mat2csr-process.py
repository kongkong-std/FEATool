#!/usr/bin/env python3
import argparse
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, issparse

def save_csr(matrix, filename):
    """保存矩阵为自定义 CSR 格式"""
    if not issparse(matrix):
        matrix = csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    nrows, ncols = matrix.shape
    nnz = matrix.nnz

    with open(filename, "w") as f:
        f.write(f"{nrows} {ncols} {nnz}\n")
        for r in matrix.indptr:
            f.write(f"{r}\n")
        for c in matrix.indices:
            f.write(f"{c}\n")
        for v in matrix.data:
            f.write(f"{v:021.16e}\n")

def save_vector(vec, filename):
    """保存向量为自定义格式"""
    vec = np.asarray(vec).squeeze()
    n = vec.shape[0]
    with open(filename, "w") as f:
        f.write(f"{n}\n")
        for v in vec:
            f.write(f"{v:021.16e}\n")

def load_mat_variable(filename, varname):
    """
    自动支持 MATLAB v7.2 及以下 (scipy) 和 v7.3 (HDF5)
    支持稀疏矩阵和普通矩阵/向量
    """
    try:
        # 尝试 v7.2 及以下
        from scipy.io import loadmat
        data = loadmat(filename)
        if varname not in data:
            raise KeyError(f"Variable {varname} not found in {filename}")
        return data[varname]

    except NotImplementedError:
        # v7.3 HDF5 格式
        import h5py
        with h5py.File(filename, "r") as f:
            if varname not in f:
                raise KeyError(f"Variable {varname} not found in {filename}")
            obj = f[varname]

            if isinstance(obj, h5py.Group):
                # 检查是否为稀疏矩阵
                if all(k in obj for k in ("data", "ir", "jc")):
                    data = np.array(obj["data"])
                    ir = np.array(obj["ir"])
                    jc = np.array(obj["jc"])
                    # 推算 nrows, ncols
                    ncols = len(jc) - 1
                    nrows = ir.max() + 1
                    mat = csc_matrix((data, ir, jc), shape=(nrows, ncols))
                    return mat
                else:
                    # 如果不是稀疏矩阵 group，尝试取第一层 dataset
                    raise ValueError(f"{varname} is a group but not a recognized sparse matrix")

            elif isinstance(obj, h5py.Dataset):
                # 普通矩阵/向量
                return np.array(obj[()])
            else:
                raise TypeError(f"Unsupported HDF5 object type for {varname}: {type(obj)}")

def main():
    parser = argparse.ArgumentParser(description="Convert MATLAB .mat matrix/vector to custom CSR or vector format.")
    parser.add_argument("--mat", required=True, help="path to .mat file")
    parser.add_argument("--var", required=True, help="variable name of .mat data")
    parser.add_argument("--out", required=True, help="path to output file")
    parser.add_argument("--type", required=True, choices=["matrix", "vector"], help="data type: matrix or vector")
    args = parser.parse_args()

    arr = load_mat_variable(args.mat, args.var)

    if args.type == "matrix":
        save_csr(arr, args.out)
    elif args.type == "vector":
        save_vector(arr, args.out)

if __name__ == "__main__":
    main()

"""
用法示例：
python mat2csr-process.py --mat matrix.mat --var A --out a.txt --type matrix
python mat2csr-process.py --mat vector.mat --var b --out b.txt --type vector
"""
