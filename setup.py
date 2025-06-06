from setuptools import setup, find_packages

setup(
    name='CSA_S16_python',
    version='0.1.1',
    packages=find_packages(),
    install_requires=[
        'handcalcs',
        'forallpeople',
        'IPython',
    ],
    author='Sam Richardson',
    author_email='sam.richardson@mai.utoronto.ca',
    description='Steel design equations based on CSA S16.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/samrich98/CSA_S16_python.git',  # Change this if you publish it
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.7',
)
