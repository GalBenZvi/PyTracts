import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyTracts",  # Replace with your own username
    version="0.0.1",
    author="Gal Hershkovitz",
    author_email="hershkovitz1@mail.tau.ac.il",
    description="Automated tractography pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Hershkovitz-hub/PyTracts",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    data_files=["logs", "atlases"],
)
