from setuptools import setup, find_packages

with open("wisp_light/requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="wisp_light",
    version="0.1.0",
    author="Maxwell",
    description="Un outil pour le traitement et l'analyse de données génomiques avec apprentissage automatique",
    long_description=open("wisp_light/README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["notebooks", "tests", "__pycache__"]),

    install_requires=requirements,
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    entry_points={
        "console_scripts": [
            "tk-infer = wisp_light.prediction.tk_infer:main",
        ]
    },
)