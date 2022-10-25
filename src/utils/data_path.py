import os

# Taking DATA_PATH via command line arguments and then passing it to all
# modules is tricky, so instead I pass it via environment variables
if "DATA_PATH" in os.environ:
    DATA_PATH = os.environ["DATA_PATH"]
else:
    _SCRIPT_DIR = os.path.dirname(__file__)
    DATA_PATH = os.path.abspath(os.path.join(_SCRIPT_DIR, "../../data/"))
    # Canary in case I move files around and forget to change the path
    assert os.path.isdir(DATA_PATH)
    assert os.path.isdir(DATA_PATH + "/original")
    # assert os.path.isdir(DATA_PATH + "/interim")
    assert os.path.isdir(DATA_PATH + "/database")

def get_path(s: str, create=False, exist_ok=True) -> str:
    path = os.path.join(DATA_PATH, s)
    if create:
        os.makedirs(path, exist_ok = exist_ok)
    else:
        assert(os.path.isdir(path)), f"Path {path} does not exist"
    return path


