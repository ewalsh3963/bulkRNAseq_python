import os
import glob

def ensure_directory(path, critical=False):
    if not os.path.isdir(path):
        if critical:
            raise FileNotFoundError(f"Critical directory does not exist: {path}")
        else:
            os.makedirs(path)

def find_files(parent_dir, extension, wd=None, recursive=False):
    """
    Searches for files within a specified directory using glob patterns.

    Args:
        parent_directory (str): The starting directory for the search.
                                This can be an absolute or relative path.
        extension (str): The file extension to search for (e.g., '.txt', '.py', 'jpg').
                         Include the dot if you want to match files with an extension.
                         Use an empty string '' or '*' to match any extension (if combined with *).
        working_directory (str, optional): An optional working directory to change into
                                           before performing the search. If None, the
                                           current working directory is used.
        recursive (bool): If True, the search will include subdirectories recursively.
                          If False, only the parent_directory is searched.

    Returns:
        list: A list of paths to the found files, relative to the working_directory
              or the directory where the script is run if no working_directory is specified.
    """
    original_cwd = os.getcwd()  # Store the original current working directory
    found_files = []

    try:
        if wd:
            os.chdir(wd)
            print(f"Changed working directory to: {os.getcwd()}")
        else:
            print(f"Searching from current working directory: {os.getcwd()}")

        # Ensure the extension starts with a dot if it's not empty and doesn't already have one
        if extension and not extension.startswith('.'):
            search_extension = f"*.{extension}"
        elif extension == "": # If extension is empty, match any file name (e.g., "file." for files with no extension)
            search_extension = "*" # Match any file
        else: # e.g. '.txt' or '*'
            search_extension = f"*{extension}"

        # Construct the glob pattern
        if recursive:
            # Use '**' to match any directories and subdirectories recursively
            # os.path.join handles path separators correctly across OS
            search_pattern = os.path.join(parent_dir, '**', search_extension)
            found_files = glob.glob(search_pattern, recursive=True)
        else:
            # Search only in the parent_directory
            search_pattern = os.path.join(parent_dir, search_extension)
            found_files = glob.glob(search_pattern, recursive=False)

    except FileNotFoundError:
        print(f"Error: Directory '{wd or os.getcwd()}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Always change back to the original working directory
        if os.getcwd() != original_cwd:
            os.chdir(original_cwd)
            print(f"Changed back to original working directory: {os.getcwd()}")

    return found_files