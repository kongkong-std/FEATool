#!/usr/bin/env python3
import argparse
import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix, issparse

def save_csr(matrix, filename):
    
    if not issparse(matrix):
        matrix = csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    nrows, ncols = matrix.shape
    nnz = matrix.nnz

    with open(filename, "w") as f:
        # header
        f.write(f"{nrows} {ncols} {nnz}\n")
        
        # row_ptr
        for r in matrix.indptr:
            f.write(f"{r}\n")

        # col_idx
        for c in matrix.indices:
            f.write(f"{c}\n")

        # val
        for v in matrix.data:
            f.write(f"{v:021.16e}\n")


def save_vector(vec, filename):
    
    vec = np.asarray(vec).squeeze()
    n = vec.shape[0]

    with open(filename, "w") as f:
        # header
        f.write(f"{n}\n")
        # values
        for v in vec:
            f.write(f"{v:021.16e}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert MATLAB .mat matrix/vector to custom CSR or vector format.")
    parser.add_argument("--mat", required=True, help="path to .mat file")
    parser.add_argument("--var", required=True, help="variable name of .mat data")
    parser.add_argument("--out", required=True, help="path to output")
    parser.add_argument("--type", required=True, choices=["matrix", "vector"], help="data type: matrix or vector")
    args = parser.parse_args()

    data = loadmat(args.mat)
    if args.var not in data:
        raise KeyError(f"variable {args.var} does not exist {args.mat}")

    arr = data[args.var]

    if args.type == "matrix":
        save_csr(arr, args.out)
    elif args.type == "vector":
        save_vector(arr, args.out)


if __name__ == "__main__":
    main()

'''
python mat2csr.py --mat matrix.mat --var A --out a.txt --type matrix
python mat2csr.py --mat vector.mat --var b --out b.txt --type vector
'''
