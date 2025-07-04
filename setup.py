from setuptools import setup, find_packages


with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='fluhp',
    version='1.0.0',
    author='viralInformatics',
    author_email='pys2013@hnu.edu.cn',
    description="A bioinformatics tool for host prediction of emerging influenza A viruses using genomic sequences.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/viralInformatics/FluHostPredictor',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'fluhp=fluhp.fluhp:main',
        ],
    },
    package_data={
        '': ['data/*', 'model/*','test/*', '18Mid/*', 'temp/*','app/*', 'querySeq/*','result/*'],
    },
    install_requires=required,
    zip_safe=False,
 
    classifiers=[
        'Programming Language :: Python :: 3.6',
        # 'Programming Language :: Python :: 3.7',
        # 'Programming Language :: Python :: 3.8',
        # 'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
    ],
    # python_requires='>=3.6',
    python_requires='>=3.6, <3.8',
    keywords = 'influenza virus, risk assessment,genome'
)
