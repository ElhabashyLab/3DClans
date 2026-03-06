import logging
import os
import shutil
import requests

logger = logging.getLogger(__name__)


def download_file(url, output_path):
    """
    Downloads a file from a given URL and saves it to a specified local path.
    """
    try:
        response = requests.get(url, stream=True, timeout=10)
        # raise exception if request failed
        response.raise_for_status()
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'wb') as f:
            # download file in 8KB chunks
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        logger.warning(f"Failed to download {url}: {e}")
        return False
    

def reset_dir_content(dir_path):
    """
    Deletes the content of the specified directory and creates it if it does not exists.
    """
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    
    
def copy_dir_content(source_dir, target_dir):
    """
    Copies all files and subdirectories from source_dir into target_dir.
    Creates target_dir if it does not exist.

    Args:
        source_dir (str): Path to the directory to copy FROM
        target_dir (str): Path to the directory to copy TO
    """
    # Ensure target exists
    os.makedirs(target_dir, exist_ok=True)

    # Iterate through items in source_dir
    for item in os.listdir(source_dir):
        src = os.path.join(source_dir, item)
        dst = os.path.join(target_dir, item)

        if os.path.isdir(src):
            # Copy a directory recursively
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            # Copy a single file
            shutil.copy2(src, dst)
