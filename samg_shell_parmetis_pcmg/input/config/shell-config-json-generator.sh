#!/bin/bash

SCRIPT="config-json-generator.py"

# -----------------------------------------------------
# 配置区（你只需要维护这里）
# -----------------------------------------------------

# 全局参数（可维护）
BASE="/home/kongkong/Documents/DocFEMTool/FEATool/samg_shell_parmetis_pcmg"
NP="np-4"
PRE_SMOOTH=2
POST_SMOOTH=2
NUM_LEVEL=10
NUM_COARSE_VTX=100
EST_SIZE_AGG=4
PS_NUM_STEPS=1
PS_TYPE=0
PS_SCALE=0.67

# 模型名数组
MODELS=("flat-shell" "inclined-shell" "spherical-shell")

# 每个模型的 DOF（保持一一对应）
declare -A DOF=(
    ["flat-shell"]=480
    ["inclined-shell"]=300
    ["spherical-shell"]=684
)

# -----------------------------------------------------
# 执行区（无需修改）
# -----------------------------------------------------

for model in "${MODELS[@]}"; do
    echo "============================================"
    echo "  Processing model: ${model}"
    echo "============================================"

    model_dof=${DOF[$model]}

    INPUT_DIR="${BASE}/input/${model}/dof-${model_dof}/${NP}"
    OUTPUT_JSON="${BASE}/input/config/comsol-model/${model}/dof-${model_dof}.json"

    python "$SCRIPT" \
        --output_json "$OUTPUT_JSON" \
        --file_mat "$INPUT_DIR" \
        --file_rhs "$INPUT_DIR" \
        --file_vtx "$INPUT_DIR" \
        --file_adj "$INPUT_DIR" \
        --pre_smooth $PRE_SMOOTH \
        --post_smooth $POST_SMOOTH \
        --num_level $NUM_LEVEL \
        --num_coarse_vtx $NUM_COARSE_VTX \
        --est_size_agg $EST_SIZE_AGG \
        --ps_num_steps $PS_NUM_STEPS \
        --ps_type $PS_TYPE \
        --ps_scale $PS_SCALE

    echo ""
done

echo ">>> All config-json generation tasks finished."
