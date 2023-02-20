#!/usr/bin/env python
import git
import sys
import os
import yaml
from random import randint

MAX_SEED = 4_294_967_295 # u32::MAX


def bail(msg, error):
    print(f"[ERROR]: {msg}", file=sys.stderr)
    raise RuntimeError(error)


def get_git_root(path):
    from git.exc import InvalidGitRepositoryError
    """
    Attempt to obtain the Hash of a directory's git repository.
    Errors if 'path' is not a valid git directory.
    """

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        commit   = git_repo.git.rev_parse("HEAD")
        name     = os.path.basename(git_repo.git.rev_parse("--show-toplevel"))
    except InvalidGitRepositoryError as e:
        bail(f"'{path}' does not appear to be a valid git repository.", e)
    return {'name': name, 'commit': commit}


def dump_yaml(yaml_dictionary, output_path):
    with open(output_path, 'w') as file:
        documents = yaml.dump(yaml_dictionary, file)


def parse_git_commits(yaml_dictionary):
    gitted_directories = {
        "aDNA-pipeline": ".",
        "GRUPS-rs": "workflow/scripts/grups/"
    }

    for (repo, path) in gitted_directories.items():
        repo = get_git_root(path)
        yaml_dictionary["git-versions"][repo['name']] = repo['commit']


def main():
    # ---- Prepare a yaml dictionnary
    yaml_dict = {"git-versions": {} }

    # ---- Keep track of the git commit for these directories
    parse_git_commits(yaml_dict)

    # ---- Create a random "backup" seed. This can in turn be used
    #      by any tool for which the user did not explicitly request
    #      a seed.
    yaml_dict["seed"] = randint(0, MAX_SEED)

    # ---- Write all contents in the output
    dump_yaml(yaml_dict, snakemake.output.metadata)


if __name__ == '__main__':
    main()
