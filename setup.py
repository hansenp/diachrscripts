import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="diachromatic", # Replace with your own username
    version="0.0.1",
    author="Peter Hansen, Peter N Robinson",
    author_email="peter.hansen@jax.org,peter.robinson@jax.org",
    description="A package for working with capture Hi-C/Hi-C paired-end data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TheJacksonLaboratory/diachrscripts",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
    python_requires='>=3.8',
)