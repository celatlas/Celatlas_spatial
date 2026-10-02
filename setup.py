import setuptools
import os
import re

# Read version from __init__.py without importing
def get_version():
    init_path = os.path.join(os.path.dirname(__file__), 'celatlas_spatial', '__init__.py')
    with open(init_path, 'r', encoding='utf-8') as f:
        content = f.read()
        version_match = re.search(r'^__VERSION__\s*=\s*["\']([^"\']*)["\']', content, re.MULTILINE)
        if version_match:
            return version_match.group(1)
        raise RuntimeError("Unable to find version string.")

def get_assay_list():
    init_path = os.path.join(os.path.dirname(__file__), 'celatlas_spatial', '__init__.py')
    with open(init_path, 'r', encoding='utf-8') as f:
        content = f.read()
        # Match ASSAY_LIST = ["rna", ...]
        assay_match = re.search(r'^ASSAY_LIST\s*=\s*\[(.*?)\]', content, re.MULTILINE | re.DOTALL)
        if assay_match:
            assay_str = assay_match.group(1)
            # Extract quoted strings
            assays = re.findall(r'["\']([^"\']+)["\']', assay_str)
            return assays
        return ["rna"]  # Default fallback

__VERSION__ = get_version()
ASSAY_LIST = get_assay_list()

# Read requirements.txt
with open('requirements.txt', encoding='utf-8') as fp:
    install_requires = fp.read().strip().split('\n')

entrys = [
    'celatlas=celatlas_spatial.runner.public_cli:main',
    'celatlas_spatial=celatlas_spatial.celatlas:main',
    'celatlas_spatial_batch=celatlas_spatial.runner.cli:batch_main',
    'celatlas_spatial_runner=celatlas_spatial.runner.cli:main',
    'celatlas_spatial_report=celatlas_spatial.tools.spatial_report_generator:main',
    'celatlas_scrna_report=celatlas_spatial.tools.scrna_report_generator:main',
]
for assay in ASSAY_LIST:
    entrys.append(f'multi_{assay}=celatlas_spatial.{assay}.multi_{assay}:main')
entry_dict = {
    'console_scripts': entrys,
}

setuptools.setup(
    name='celatlas-spatial',
    version=__VERSION__,
    author='Qingsong Li',
    author_email='lqs60667106@gmail.com',
    description='Spatial Transcriptomics Analysis Pipelines',
    packages=setuptools.find_packages(exclude=('tests', 'tests.*', 'private', 'private.*')),
    python_requires='>=3.11',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
    scripts=['Celatlas.sh', 'Celatlas_FFPE.sh', 'Celatlas_SN.sh', 'Celatlas_reanalysis.sh', 'Celatlas_run.sh'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
