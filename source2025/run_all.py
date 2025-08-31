import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_script(script_path: Path, extra_args=None, cwd: Path = None, add_to_pythonpath: list[str] | None = None):
    base = Path(__file__).parent
    try:
        rel = script_path.relative_to(base)
    except Exception:
        rel = script_path
    print(f"\n=== Running {rel} ===")
    if cwd:
        print(f"[INFO] Working directory: {cwd}")
    env = os.environ.copy()
    if add_to_pythonpath:
        sep = os.pathsep
        env["PYTHONPATH"] = sep.join(add_to_pythonpath + [env.get("PYTHONPATH", "")]).strip(sep)
        print(f"[INFO] PYTHONPATH: {env['PYTHONPATH']}")
    cmd = [sys.executable, str(script_path)]
    if extra_args:
        cmd.extend(extra_args)
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None, env=env)


def parse_multi_csv(s: str):
    return [x.strip().lower() for x in s.split(",") if x.strip()]


def main():
    parser = argparse.ArgumentParser(description="Run Q1–Q4 and optional visualizations.")
    parser.add_argument(
        "--visualize", "--viz",
        help="Comma-separated list: q3, q4, all (e.g., --visualize q3,q4)",
    )
    parser.add_argument(
        "--stages",
        choices=["before", "after", "both"],
        default="both",
        help="For Q3 visualization: which stages to render (default: both)",
    )
    args, unknown = parser.parse_known_args()

    base = Path(__file__).parent
    viz_dir = base / "visualize"

    for name in ["q1.py", "q2.py", "q3.py", "q4.py"]:
        p = base / name
        if p.exists():
            try:
                run_script(p, extra_args=unknown, cwd=base, add_to_pythonpath=[str(base)])
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] {name} exited with {e.returncode}")
                sys.exit(e.returncode)
        else:
            print(f"[SKIP] {name} not found at {p}")

    if not args.visualize:
        return

    requested = set()
    for tok in parse_multi_csv(args.visualize):
        if tok == "all":
            requested.update(["q3", "q4"])
        elif tok in ("q3", "q4"):
            requested.add(tok)
        else:
            print(f"[WARN] Unknown visualize target: {tok}")

    print(f"\n=== Visualization targets: {', '.join(sorted(requested)) or '(none)'} ===")

    if "q3" in requested:
        q3_viz = viz_dir / "visualize_q3.py"
        if not q3_viz.exists():
            print(f"[MISS] Q3 viz not found: {q3_viz}")
        else:
            stages = ["before", "after"] if args.stages == "both" else [args.stages]
            print(f"[OK] Q3 stages: {', '.join(stages)}")
            for stage in stages:
                try:
                    run_script(
                        q3_viz,
                        extra_args=[stage] + unknown,
                        cwd=base,
                        add_to_pythonpath=[str(base), str(viz_dir)],
                    )
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] visualize_q3.py ({stage}) exited with {e.returncode}")
                    sys.exit(e.returncode)

    if "q4" in requested:
        q4_viz = viz_dir / "visualize_q4.py"
        if not q4_viz.exists():
            print(f"[MISS] Q4 viz not found: {q4_viz}")
        else:
            print("[OK] Running Q4 visualization with NO arguments (cwd=project root)")
            try:
                run_script(
                    q4_viz,
                    extra_args=None,
                    cwd=base,
                    add_to_pythonpath=[str(base), str(viz_dir)],
                )
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] visualize_q4.py exited with {e.returncode}")
                sys.exit(e.returncode)


if __name__ == "__main__":
    main()
