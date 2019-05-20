from setuptools import setup

setup(name='wavecalc',
      version='0.1',
      description='A tool for EM wave interactions at interfaces',
      long_description=open('README.txt').read(),
      install_requires=['numpy'],
      url='https://github.com/rmgoetz/wavecalc',
      author='Ryan Goetz',
      author_email='ryan.m.goetz@gmail.com',
      license='MIT',
      packages=['wavecalc'],
      zip_safe=False)
