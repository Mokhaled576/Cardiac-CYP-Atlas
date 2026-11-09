"""Rasterize the rendered composites and panels as 600-dpi RGB PNGs."""
import argparse
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PIL import Image


def rasterize(pdf, executable):
    png = pdf.with_suffix(".png")
    subprocess.run(
        [executable, "-r", "600", "-png", "-singlefile", str(pdf), str(png.with_suffix(""))],
        check=True,
    )
    with Image.open(png) as image:
        if image.mode != "RGB" or any(abs(dpi - 600) >= 0.1 for dpi in image.info.get("dpi", (0, 0))):
            raise ValueError(f"Unexpected image mode or resolution for {png}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Root of the rendered figure folders")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    executable = shutil.which("pdftoppm")
    if not executable:
        parser.error("Install Poppler and put pdftoppm on PATH")
    folders = ("Main_Figures", "Extended_Data/Figures", "Individual_Panels")
    files = sorted(pdf for folder in folders for pdf in (args.output / folder).rglob("*.pdf"))
    if not files:
        parser.error("No rendered PDFs found")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(lambda pdf: rasterize(pdf, executable), files))
    print(f"Verified {len(files)} RGB PNGs at 600 dpi.")


if __name__ == "__main__":
    main()
