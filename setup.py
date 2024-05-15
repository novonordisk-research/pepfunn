"""Setup pepfunn."""
from __future__ import absolute_import, print_function

import io
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_namespace_packages, setup


def _read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ) as fh:
        return fh.read()


def get_version():
    with io.open(
        join(dirname(__file__), "src/pepfunn/__init__.py"), encoding="utf-8"
    ) as vf:
        for line in vf.readlines():
            if "__version__" in line:
                return line.split("=")[1].strip().replace('"', "")


with open("requirements.txt", "r") as fp:
    reqs = fp.read().split("\n")

setup(
    name="pepfunn",
    version=get_version(),
    url="",
    author="Molecular AI - Novo Nordisk",
    author_email="RAOC@novonordisk.com",
    description="Protocols for the analysis of peptides using cheminformatics and bioinformatics tools",
    long_description="%s\n%s"
    % (
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
            "", _read("README.md")
        ),
        ""
    ),
    long_description_content_type="text/markdown",
    packages=find_namespace_packages("src"),
    package_dir={"": "src"},
    package_data={'pepfunn': ['data/*']},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=reqs,
)
