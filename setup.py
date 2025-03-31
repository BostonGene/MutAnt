import setuptools

with open("requirements.txt", "r") as f:
    install_requires = f.read().split("\n")

__version__ = "0.0.0"
exec(open("mutant/__init__.py").read())

setuptools.setup(
    name="mutant",
    version=__version__,
    packages=setuptools.find_packages(),
    url="",
    license="",
    author="BostonGene",
    author_email="",
    description="",
    setup_requires=["cython", "cmake", "wheel", "setuptools"],
    install_requires=install_requires,
    include_package_data=True,
    dependency_links=["https://nexus.devbg.us/repository/pypi-all/simple"],
)
