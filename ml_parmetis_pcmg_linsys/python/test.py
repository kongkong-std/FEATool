import numpy as np
from scipy.sparse import csr_matrix

filename = "input/KKK.mat"

with open(filename, "rb") as f:
    # 读取矩阵尺寸和非零元数量（假设 8 bytes integer）
    num_row = np.fromfile(f, dtype=np.int64, count=1)[0]
    num_col = np.fromfile(f, dtype=np.int64, count=1)[0]
    num_nonzeros = np.fromfile(f, dtype=np.int64, count=1)[0]

    print(f"num_row={num_row}, num_col={num_col}, num_nonzeros={num_nonzeros}")

    # 读取 indptr，长度 = num_row + 1, 4 bytes integer
    indptr = np.fromfile(f, dtype=np.int32, count=num_row + 1)
    print("indptr[:10] =", indptr[:10])

    # 读取 indices，长度 = num_nonzeros, 4 bytes integer
    indices = np.fromfile(f, dtype=np.int32, count=num_nonzeros)
    print("indices[:10] =", indices[:10])

    # 读取 data，长度 = num_nonzeros, 8 bytes double
    data = np.fromfile(f, dtype=np.float64, count=num_nonzeros)
    print("data[:10] =", data[:10])

# 构建 CSR 矩阵
A = csr_matrix((data, indices, indptr), shape=(num_row, num_col))
print("CSR matrix shape:", A.shape)
print("CSR matrix nnz:", A.nnz)

# 可以打印前几行稀疏矩阵内容
for i in range(min(10, num_row)):
    row_data = A[i].toarray().flatten()
    print(f"row {i}:", row_data)
