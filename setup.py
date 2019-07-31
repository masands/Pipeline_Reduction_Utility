import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="example-pkg-your-username",
    version="0.0.1a1",
    author="Sahil Masand",
    author_email="masand.sahil@gmail.com",
    description="A toolkit to for Flow Assurance engineers to reduce the number of sections in the flowline bathymetry for increased simulation speed, whilst maintaining important angle classes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)