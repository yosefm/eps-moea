from distutils.core import setup

setup(
    name='eps_moea',
    version='1.0b',
    description='Python implementation of the Epsilon MOEA algorithm',
    author='Yosef Meller',
    author_email='mellerf@netvision.net.il',
    url='http://wiki.github.com/yosefm/eps-moea',
    package_dir = {'eps_moea': 'py_eps_moea'},
    packages=['eps_moea'], 
)
