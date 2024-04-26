from setuptools import setup, find_packages
  
with open('requirements.txt') as f:
    requirements = f.readlines()
  
long_description = 'Compare ATAC-seq datasets against loci datasets to see if there is a specific enrichment between a particular loci set in a cell type'
  
setup(
        name ='gaiaAssociation',
        version ='1.2.4',
        author ='Samuel Rosean',
        author_email ='samuel.rosean@einsteinmed.edu',
        url ='https://github.com/samrosean/gaiaAssociation',
        description ='Compare ATAC-seq data to loci.',
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='MIT',
        packages = find_packages(),
        entry_points ={
            'console_scripts': [
                'gaia = gaiaAssociation.gaiaAssociation:main'
            ]
        },
        classifiers =[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        keywords ='ATAC-seq loci enrichment',
        install_requires = requirements,
        zip_safe = False
)
