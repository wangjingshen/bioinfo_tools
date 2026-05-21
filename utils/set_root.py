import sys
from pathlib import Path

def add_root(marker="utils"):
    """自动向上找包含 marker 文件夹的根目录，并加入 sys.path"""
    root = Path(__file__).resolve().parent
    while not (root / marker).exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError(f"{marker} 文件夹未找到！")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
pipeline_root = Path(__file__).resolve().parents[1]