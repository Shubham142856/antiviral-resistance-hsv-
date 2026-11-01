import sys

modules = ['requests', 'Bio', 'scipy', 'sklearn', 'torch', 'transformers']
for mod in modules:
    try:
        m = __import__(mod)
        ver = getattr(m, '__version__', 'unknown')
        print(f"{mod}: available (version {ver})")
    except ImportError as e:
        print(f"{mod}: NOT available ({e})")

