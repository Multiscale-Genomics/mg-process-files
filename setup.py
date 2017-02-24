from setuptools import setup

setup(
    name='mg-process-tsv',
    packages=['mg-process-tsv'],
    include_package_data=True,
    install_requires=[
        'numpy', 'h5py', 'pycompss'
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
)
