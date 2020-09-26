from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# long_description(後述)に、GitHub用のREADME.mdを指定
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'pybeads',
    version = '1.0',
    license = 'MIT',
    author = 'Kotaro Saito',
    author_email = 'kotaro.saito@kek.jp',
    description = 'Baseline estimation and denoizing algorithm using sparcity',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/skotaro/pybeads',
    packages = find_packages(exclude=['sample_data']),
    classifiers = [
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        "Operating System :: OS Independent",
    ], # ref. https://pypi.org/classifiers/
    python_requires='>=3.6',
    install_requires = ['numpy', 'scipy'],
    keywords = 'beads baseline background',
)