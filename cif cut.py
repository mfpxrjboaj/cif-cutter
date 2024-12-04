from ase.io import read, write
from ase import Atoms
import numpy as np
import os

# 加载结构
cif_file = r"E:\project pp\练习用\cif分割\11.cif"
structure = read(cif_file)

# 定义分割数量 nx, ny, nz，确保 nx * ny * nz = 10
nx, ny, nz = 5, 2, 1
if nx * ny * nz != 10:
    raise ValueError("nx, ny, nz 的乘积必须为 10！")

# 获取原始晶胞信息
cell = structure.cell  # 晶胞矩阵
positions = structure.positions  # 原子坐标
scaled_positions = structure.get_scaled_positions()  # 原子归一化坐标

# 创建输出文件夹
output_folder = r"E:\project pp\练习用\cif分割\split_structures"
os.makedirs(output_folder, exist_ok=True)

count = 0
# 分割晶胞
for ix in range(nx):
    for iy in range(ny):
        for iz in range(nz):
            # 定义当前分割晶胞的范围
            min_frac = np.array([ix / nx, iy / ny, iz / nz])
            max_frac = np.array([(ix + 1) / nx, (iy + 1) / ny, (iz + 1) / nz])

            # 筛选属于当前分割区域的原子
            mask = np.all((scaled_positions >= min_frac) & (scaled_positions < max_frac), axis=1)
            sub_positions = positions[mask]
            sub_symbols = np.array(structure.get_chemical_symbols())[mask]

            # 创建子晶胞
            subcell = cell.copy()
            subcell[0] /= nx
            subcell[1] /= ny
            subcell[2] /= nz
            sub_structure = Atoms(
                symbols=sub_symbols,
                positions=sub_positions,
                cell=subcell,
                pbc=True
            )

            # 保存子晶胞
            output_file = os.path.join(output_folder, f"split_{count}.cif")
            write(output_file, sub_structure)
            count += 1

print(f"分割完成！所有子结构已保存到: {output_folder}")

