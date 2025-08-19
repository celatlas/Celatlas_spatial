import setuptools
from celatlas_spatial.__init__ import __VERSION__, ASSAY_LIST


with open('requirements.txt') as fp:
    install_requires = fp.read()

entrys = ['celatlas_spatial=celatlas_spatial.celatlas:main',]
for assay in ASSAY_LIST:
    entrys.append(f'multi_{assay}=celatlas_spatial.{assay}.multi_{assay}:main')
entry_dict = {
    'console_scripts': entrys,
}

setuptools.setup(
    name='celatlas-spatial',
    version=__VERSION__,
    author='Celatlas',
    author_email='rd@celatlas.com; lqs60667106@gmail.com',
    description='Spatial Transcriptomics Analysis Pipelines',
    packages=setuptools.find_packages(),
    python_requires='=3.9',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
    scripts=['Celatlas.sh'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
