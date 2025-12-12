#!/bin/bash

PROCESS_SCRIPT="mesh-process.py"
INPUT_DIR="input"

# 声明关联数组：模型名 → mesh文件名
declare -A MESH_FILE=(
    ["flat-shell"]="mesh-flat-shell.mphtxt"
    ["inclined-shell"]="mesh-inclined-shell.mphtxt"
    ["spherical-shell"]="mesh-curved-shell.mphtxt"
)

MODELS=("flat-shell" "inclined-shell" "spherical-shell")

for model in "${MODELS[@]}"; do
    echo "============================================"
    echo "  Processing: ${model}"
    echo "============================================"

    mesh_name="${MESH_FILE[$model]}"

    python "$PROCESS_SCRIPT" \
        --mesh          "${INPUT_DIR}/${model}/${mesh_name}" \
        --vertex_coor   "${INPUT_DIR}/${model}/vertex_coor.txt" \
        --adjacent_list "${INPUT_DIR}/${model}/adjacent_list.txt" \
        --entity        "${INPUT_DIR}/${model}/entity.txt"

    echo ""
done

echo ">>> All tasks completed."
