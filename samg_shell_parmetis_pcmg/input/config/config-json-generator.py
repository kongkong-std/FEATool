import argparse
import json
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="config json file generator"
    )

    # -------- file paths --------
    parser.add_argument("--file_mat", type=str, required=True,
                        help="Path to distributed matrix file")
    parser.add_argument("--file_rhs", type=str, required=True,
                        help="Path to distributed RHS file")
    parser.add_argument("--file_vtx", type=str, required=True,
                        help="Path to distributed vertex coordinate file")
    parser.add_argument("--file_adj", type=str, required=True,
                        help="Path to distributed adjacency list file")

    # -------- mg parameters --------
    parser.add_argument("--pre_smooth", type=int, default=2,
                        help="Number of pre-smoothing iterations")
    parser.add_argument("--post_smooth", type=int, default=2,
                        help="Number of post-smoothing iterations")
    parser.add_argument("--num_level", type=int, default=10,
                        help="Number of multigrid levels")
    parser.add_argument("--num_coarse_vtx", type=int, default=100,
                        help="Number of coarse-level vertices")
    parser.add_argument("--est_size_agg", type=int, default=4,
                        help="Estimation size of each aggregate")
    parser.add_argument("--ps_num_steps", type=int, default=1,
                        help="Number of prolongation operator smoothed steps")
    parser.add_argument("--ps_type", type=int, default=0,
                        help="Prolongation operator smoother type")
    parser.add_argument("--ps_scale", type=float, default=0.67,
                        help="Scaling prolongation smoother")

    # -------- output json path --------
    parser.add_argument("--output_json", type=str, required=True,
                        help="Output path for generated config json file")

    return parser.parse_args()


def generate_json(args):
    config = {
        "file": {
            "file_mat": args.file_mat,
            "file_rhs": args.file_rhs,
            "file_vtx": args.file_vtx,
            "file_adj": args.file_adj
        },
        "mg": {
            "pre_smooth": args.pre_smooth,
            "post_smooth": args.post_smooth,
            "num_level": args.num_level,
            "num_coarse_vtx": args.num_coarse_vtx,
            "est_size_agg": args.est_size_agg,
            "ps_num_steps": args.ps_num_steps,
            "ps_type": args.ps_type,
            "ps_scale": args.ps_scale
        }
    }

    with open(args.output_json, "w") as f:
        json.dump(config, f, indent=4)

    print(f"[OK] JSON config written to: {args.output_json}")


if __name__ == "__main__":
    args = parse_args()
    generate_json(args)


# usage
# python config-json-generator.py \
#     --output_json </output/path/for/json/file> \
#     --file_mat </base/path/to/distributed/mat/file> \
#     --file_rhs </base/path/to/distributed/rhs/file> \
#     --file_vtx </base/path/to/distributed/vtx/file> \
#     --file_adj </base/path/to/distributed/adj/file> \
#     --pre_smooth <presmoothing times> \
#     --post_smooth <postsmoothing times> \
#     --num_level <number of levels> \
#     --num_coarse_vtx <number of coarse lelve vertices> \
#     --est_size_agg <estimation of each aggregate>
#
