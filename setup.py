from distutils.core import setup


requires = ['numpy','copy','time','scipy']

packages = [
		'mkin4py',
		'mkin4py.solver',
]

package_dir = {'mkin4py' : 'mkin4py'}
package_data = { 'mkin4py' : []}


setup(
    name='mkin4py',
    version='1.0',
    packages=packages,
    license='The MIT License (MIT)',
    author = 'Gabriel S. Gusmao',
    author_email = 'gusmaogabriels@gmail.com',
    url = 'https://github.com/gusmaogabriels/mkin4py',
    download_url = 'https://github.com/gusmaogabriels/mkin4py/tarball/v1.0',
    keywords = ['python', 'microkinetics', 'linearized', 'solver', 'chemical', 'kinetics', 'power-law'],
    package_data = package_data,
    package_dir = package_dir,

)
