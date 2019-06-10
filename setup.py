from setuptools import setup

setup(name='wavecalc',
      version='0.1',
      description='A tool for EM wave interactions at interfaces',
      long_description=open('README.txt').read(),
      install_requires=['numpy'],
      classifiers=[
            'Development Status :: 1 - Planning',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Physics'    
      ],
      url='https://github.com/rmgoetz/wavecalc',
      author='Ryan Goetz',
      author_email='ryan.m.goetz@gmail.com',
      license='MIT',
      packages=['wavecalc'],
      zip_safe=False)
