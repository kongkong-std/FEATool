#!/bin/bash

PROCESS_SCRIPT="parfile-process.py"
INPUT_DIR="input"

# -------------------------
# 可维护区：所有可配置内容（你只需要改这里）
# -------------------------

# 进程数（可随时修改）
NUM_PROC=4

# 要处理的模型列表
MODELS=(
    "flat-shell"
    "inclined-shell"
    "spherical-shell"
)

# -------------------------
# 执行区（无需修改）
# -------------------------

for model in "${MODELS[@]}"; do
    echo "============================================"
    echo "  Processing model: ${model}"
    echo "============================================"

    MODEL_DIR="${INPUT_DIR}/${model}"

    python "$PROCESS_SCRIPT" \
        --num_proc ${NUM_PROC} \
        --vertex_coor        "${MODEL_DIR}/vertex_coor.txt" \
        --output_vertex_coor "${MODEL_DIR}" \
        --adjacent_list      "${MODEL_DIR}/adjacent_list.txt" \
        --output_adjacent_list "${MODEL_DIR}" \
        --mat_file           "${MODEL_DIR}/A.txt" \
        --output_mat_file    "${MODEL_DIR}" \
        --vec_file           "${MODEL_DIR}/b.txt" \
        --output_vec_file    "${MODEL_DIR}"

    echo ""
done

echo ">>> All parfile-process tasks finished."
