import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sparce",
    version="0.0.1",
    author="Michael Skaro",
    author_email="mskaro.ms@gmail.com",
    description="A python package for automated feature selection",
    long_description_content_type="text/markdown",
    url="https://github.com/michaelSkaro/sparce/",
    project_urls={
        "Bug Tracker": "https://github.com/michaelSkaro/sparce/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #package_dir={'fs':'src/'},
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    install_requires=['seaborn','matplotlib','numpy','pandas','scikit-learn']
)
