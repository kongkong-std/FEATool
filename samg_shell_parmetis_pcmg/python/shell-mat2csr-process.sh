#!/bin/bash

PROCESS_SCRIPT="mat2csr-process.py"
INPUT_DIR="input"

# -------------------------
# 配置区：每个模型单独维护
# -------------------------

MODELS=("flat-shell" "inclined-shell" "spherical-shell")

# 矩阵 mat 文件名
declare -A MAT_FILE=(
    ["flat-shell"]="mat.mat"
    ["inclined-shell"]="mat.mat"
    ["spherical-shell"]="mat.mat"
)

# RHS 文件名
declare -A RHS_FILE=(
    ["flat-shell"]="rhs.mat"
    ["inclined-shell"]="rhs.mat"
    ["spherical-shell"]="rhs.mat"
)

# 矩阵中的变量名称（可能不是 K）
declare -A MAT_VAR=(
    ["flat-shell"]="K"
    ["inclined-shell"]="K"
    ["spherical-shell"]="K"
)

# RHS 中的变量名称（可能不是 L）
declare -A RHS_VAR=(
    ["flat-shell"]="L"
    ["inclined-shell"]="L"
    ["spherical-shell"]="L"
)

# -------------------------
# 执行区：无需修改
# -------------------------

for model in "${MODELS[@]}"; do
    echo "============================================"
    echo "  Processing model: ${model}"
    echo "============================================"

    MODEL_DIR="${INPUT_DIR}/${model}"

    # Matrix → A.txt
    python "$PROCESS_SCRIPT" \
        --mat "${MODEL_DIR}/${MAT_FILE[$model]}" \
        --var "${MAT_VAR[$model]}" \
        --out "${MODEL_DIR}/A.txt" \
        --type matrix

    # RHS → b.txt
    python "$PROCESS_SCRIPT" \
        --mat "${MODEL_DIR}/${RHS_FILE[$model]}" \
        --var "${RHS_VAR[$model]}" \
        --out "${MODEL_DIR}/b.txt" \
        --type vector

    echo ""
done

echo ">>> All mat2csr tasks completed."
