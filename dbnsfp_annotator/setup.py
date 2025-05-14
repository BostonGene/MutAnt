import setuptools

with open("requirements.txt") as f:
    install_requires = f.read().split("\n")

__version__ = "0.0.0"
exec(open("dbnsfp_annotation/__init__.py").read())

setuptools.setup(
    name="dbnsfp_annotation",
    version=__version__,
    packages=["dbnsfp_annotation"],
    url="",
    license="",
    author="BostonGene",
    author_email="",
    description="",
    setup_requires=["cython", "cmake", "wheel", "setuptools"],
    install_requires=install_requires,
    include_package_data=False,
    dependency_links=["https://nexus.devbg.us/repository/pypi-all/simple"],
)
