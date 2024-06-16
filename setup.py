import pybind11
from distutils.core import setup, Extension

ext_modules = [
    Extension(
        'cstrlbg', # название нашей либы
        ['MainPy.cpp', 'GaussianMixture.cpp', 'LindeBuzoGray.cpp'], # файлики которые компилируем
        include_dirs=[pybind11.get_include()],  # не забываем добавить инклюды pybind11
        language='c++',
        extra_compile_args=['-std=c++11'],  # используем с++11
    ),
]

setup(
    name='cstrlbg',
    version='0.1',
    author='Dmitry Lednov',
    description='Em and LBG clustering',
    ext_modules=ext_modules,
    requires=['pybind11']  # не забываем указать зависимость от pybind11
)