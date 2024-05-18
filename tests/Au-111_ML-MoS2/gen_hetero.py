# Please use the revised versions of the following packages by me to execute this script:
# https://github.com/hongyi-zhao/OgreInterface
# https://github.com/hongyi-zhao/hetbuilder/

#Problem in optimisation of H on a tmdc/metal hybrid.
#https://www.vasp.at/forum/viewtopic.php?t=19400

#https://github.com/DerekDardzinski/OgreInterface/issues/6#issuecomment-1959468808
from OgreInterface.generate import SurfaceGenerator
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot
from ase.io import read
from os import makedirs
from os.path import exists, join
from shutil import rmtree

# Define the name of the directory for storing generated files
#output_directory = "results"

# Check if the directory exists, remove it if it does, then create a new one
#if exists(output_directory):
#    rmtree(output_directory)
#makedirs(output_directory)

# 基于现有数据库的优化结果，进一步用适当算法过滤分析合理大小的超胞结构。
# 从MS中构建相应面，再看能否得到更小的匹配超胞。

# Create a surface generator instance from a CIF file
subs = SurfaceGenerator.from_file(
#    "Au_mp-81.cif",
    # 这个的111面，hetbuilder: raise Exception("Structure does not appear to be 2d.")
    "Au_mc3d_sg229.cif",
    miller_index=[1, 1, 1],  # Specify the Miller index for the surface
    layers=1,                # Number of layers in the generated surface
    vacuum=15,                # Vacuum spacing around the surface
    refine_structure=False  # Don't refine the initial structure; uncomment if needed
)


films = SurfaceGenerator.from_file(
    "MoS2_mp-1434.cif",
    miller_index=[0, 0, 1],  # Specify the Miller index for the surface
    layers=1,                # Number of layers in the generated surface
    vacuum=15,                # Vacuum spacing around the surface
    refine_structure=False  # Don't refine the initial structure; uncomment if needed
)

# Save the path for the generated file to avoid repetition
#generated_bottom_file = join(output_directory, "Au_mp-81_POSCAR_111")
#subs._slabs[0].write_file(output=generated_bottom_file, orthogonal=True)

## Read in the structure files
#bottom = read(generated_bottom_file)


#此代码段展示了如何直接从_orthogonal_slab_structure提取所需数据并构建ase.Atoms对象，从而避免了写文件的步骤。注意，此代码示例可能需要根据_orthogonal_slab_structure的实际结构进行调整。
#为了从_orthogonal_slab_structure中提取必要的信息并创建一个ase.Atoms对象，而不必将数据写入文件，您可以按照以下步骤操作：

#提取晶胞信息：从_orthogonal_slab_structure中获取晶胞参数（如abc、angles以及基矢A、B、C）。

#提取原子位置和元素符号：遍历_orthogonal_slab_structure中的PeriodicSite条目，从中提取原子的坐标位置和元素符号。

#创建ase.Atoms对象：使用提取的晶胞信息和原子位置、符号信息来构建ase.Atoms对象。

#以下是基于上述步骤的具体代码实现：

from ase import Atoms
import numpy as np

sub = subs._slabs[0]._orthogonal_slab_structure
sub_cell_matrix = sub.lattice.matrix
sub_positions = [site.coords for site in sub]
sub_symbols = [site.specie.symbol for site in sub]

film = films._slabs[0]._orthogonal_slab_structure
film_cell_matrix = film.lattice.matrix
film_positions = [site.coords for site in film]
film_symbols = [site.specie.symbol for site in film]

# See ~/.wingpro10/preferences for the related settings used by me.

# Compared to PyCharm, Wing's "Stack Data" can only provide very limited data.
# https://ask.wingware.com/question/8650/compared-to-pycharm-wings-stack-data-can-only-provide-very-limited-data/

# 创建ase.Atoms对象
bottom = Atoms(symbols=sub_symbols, positions=sub_positions, cell=sub_cell_matrix, pbc=True)
#top = Atoms(symbols=film_symbols, positions=film_positions, cell=film_cell_matrix, pbc=True)

#https://www.materialscloud.org/discover/mc2d/details/MoS2
top = read("MoS2_mc2d.cif")


# 现在可以直接使用这个atoms对象或使用上面的read方法来进一步处理：
alg = CoincidenceAlgorithm(bottom, top)

# Run the algorithm with a set of parameters
results = alg.run(tolerance=0.2)

#for i, res in enumerate(results):
#    # Save generated files in the specified directory
#    stack = join(output_directory, str(res.stack.symbols))
#    res.stack.write(stack, format='vasp')

# Set up paths for the interactive plot
iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=0.5)
iplot.plot_results()

# if True:
    # pass
    
    
    
# 创建ase.Atoms对象时，下面这些方法分别等价吗？
#In [20]: [site.coords for site in structure.sites]
#Out[20]: [array([0.      , 0.      , 8.429029])]

#In [21]: [site.coords for site in structure]
#Out[21]: [array([0.      , 0.      , 8.429029])]


#In [25]: [element.symbol for element in structure.species]
#    ...: 
#Out[25]: ['Au']

#In [26]: [site.specie.symbol for site in structure]
#Out[26]: ['Au']

#是的，基于您提供的输出，这些方法在实际使用中是等价的，但它们访问和提取所需信息的途径略有不同。让我们详细解析一下：

#提取原子坐标

#使用structure.sites: [site.coords for site in structure.sites] 这行代码通过访问structure.sites来迭代所有原子位置（或称为"站点"），并收集它们的坐标。structure.sites明确提供了一个包含所有原子站点的列表。

#直接迭代structure: [site.coords for site in structure] 这行代码显示structure对象本身支持迭代，迭代时直接返回其内部的原子站点对象，并从中提取坐标。这表明structure对象被设计为可迭代的，其迭代器直接提供对原子站点的访问。

#两种方式都有效地提取了原子的坐标信息。

#提取原子元素符号

#使用structure.species: [element.symbol for element in structure.species] 这行代码通过访问structure.species来获取一个包含所有不同原子种类的Element对象列表，然后提取这些Element对象的符号。structure.species提供了结构中每种不同元素的唯一实例。

#直接迭代structure中的原子: [site.specie.symbol for site in structure] 这行代码通过直接迭代structure对象来访问每个原子（或站点），然后提取每个原子的元素符号。这种方法逐个处理结构中的每个原子，并获取它们的元素符号。

#两种方式都成功提取了原子的元素符号，结果表明结构中包含的元素为金（Au）。虽然这些方法在操作上有所不同，但它们都达到了相同的目的。在创建ase.Atoms对象时，您可以根据实际情况选择最适合您需求的方法来提取所需的信息。

    
