#!/usr/bin/env python3
"""Rasterise 1200×630 Open Graph cards from tools/og/*-card.html via Playwright."""
from __future__ import annotations

import argparse
import socket
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CARDS = {
    "open": ("open-card.html", ROOT / "assets/og/og-open.png"),
    "research": ("research-card.html", ROOT / "assets/og/og-research.png"),
}


def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def serve(root: Path, port: int) -> subprocess.Popen[bytes]:
    return subprocess.Popen(
        [sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1"],
        cwd=root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def render(card: str) -> None:
    from playwright.sync_api import sync_playwright

    html_name, out_path = CARDS[card]
    port = free_port()
    server = serve(ROOT, port)
    try:
        time.sleep(0.4)
        url = f"http://127.0.0.1:{port}/tools/og/{html_name}"
        with sync_playwright() as p:
            browser = p.chromium.launch()
            page = browser.new_page(viewport={"width": 1200, "height": 630}, device_scale_factor=1)
            page.goto(url, wait_until="networkidle", timeout=120_000)
            page.wait_for_timeout(2500)
            page.screenshot(path=str(out_path), clip={"x": 0, "y": 0, "width": 1200, "height": 630})
            browser.close()
    finally:
        server.terminate()
        server.wait(timeout=10)
    size = __import__("struct").unpack(">II", out_path.read_bytes()[16:24])
    if size != (1200, 630):
        raise SystemExit(f"{out_path}: expected 1200×630, got {size[0]}×{size[1]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cards", nargs="*", choices=CARDS, default=list(CARDS))
    args = parser.parse_args()
    for card in args.cards:
        print(f"render {card} -> {CARDS[card][1]}", flush=True)
        render(card)


if __name__ == "__main__":
    main()
