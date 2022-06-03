from pathlib import Path
import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("build/requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="WISP",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Siegfried Dubois",
    author_email="siegfried.dubois@inria.fr",
    description="A Python application for bacterial families identification from long reads.",
    long_description=long_description,
    long_description_content_type="text/md",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=requirements,
)

# creating env path
Path(f"genomes/").mkdir(parents=True, exist_ok=True)
Path(f"genomes/train/").mkdir(parents=True, exist_ok=True)
Path(f"genomes/unk/").mkdir(parents=True, exist_ok=True)
