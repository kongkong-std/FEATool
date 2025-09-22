#!/usr/bin/env python3
import sys

def csr_to_coo(infile, outfile):
    with open(infile, 'r') as f:
        # 读取矩阵规模和非零元数
        m, n, nnz = map(int, f.readline().split())
        
        # 读取 row_pointer
        row_ptr = [int(f.readline()) for _ in range(m+1)]
        
        # 读取 col_idx
        col_idx = [int(f.readline()) for _ in range(nnz)]
        
        # 读取 val
        val = [float(f.readline()) for _ in range(nnz)]

    # 转换为 COO 格式
    rows = []
    for i in range(m):
        for j in range(row_ptr[i], row_ptr[i+1]):
            rows.append((i, col_idx[j], val[j]))

    # 写入 COO 文件
    with open(outfile, 'w') as f:
        f.write(f"{m} {n} {nnz}\n")
        for r, c, v in rows:
            f.write(f"{r} {c} {v:021.16e}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python csr2coo.py input.txt output.txt")
        sys.exit(1)
    
    csr_to_coo(sys.argv[1], sys.argv[2])
