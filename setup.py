from distutils.core import setup

setup(name='vtool',
      version='0.1',
      package_dir={'vtool': 'src'},
      packages=['vtool'],
      scripts = ["bin/g2v.py", "bin/v2g.py"]
      )
