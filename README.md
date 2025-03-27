# PDB链检测工具

这个工具用于检测PDB文件中是否存在多条链（二聚体或多聚体），特别适用于需要区分单体和二聚体结构的场景。

## 功能特点

- 检测PDB文件中的链数量
- 支持批量处理整个目录的PDB文件
- 可选递归搜索子目录
- 统计每条链的残基数量
- 生成详细的统计报告
- 支持导出结果到CSV文件
- 友好的命令行界面

## 安装依赖

```bash
pip install -r requirements.txt
```

## 使用方法

### 基本用法

```bash
python checkmonomer.py /path/to/pdb/files
```

### 带选项的用法

```bash
# 递归搜索子目录
python checkmonomer.py /path/to/pdb/files --recursive

# 导出结果到CSV文件
python checkmonomer.py /path/to/pdb/files --output results.csv

# 仅显示简要统计信息
python checkmonomer.py /path/to/pdb/files --simple
```

### 所有选项

```
usage: checkmonomer.py [-h] [--recursive] [--output OUTPUT] [--simple] directory

检测PDB文件中是否存在多条链（二聚体）

positional arguments:
  directory            PDB文件所在目录

options:
  -h, --help           show this help message and exit
  --recursive, -r      是否递归搜索子目录
  --output OUTPUT, -o OUTPUT
                       导出结果到CSV文件
  --simple, -s         仅显示简要统计信息
```

## 输出示例

脚本会生成如下输出：

1. **汇总统计**：
   - 总文件数
   - 二聚体数量
   - 单体数量
   - 处理出错数
   - 二聚体比例
   - 链数分布

2. **详细信息表格**：
   - 文件名
   - 是否二聚体
   - 链数量
   - 链ID
   - 错误信息（如有）

3. **二聚体文件列表**：
   列出所有检测到的二聚体文件名称

## CSV输出字段

导出的CSV文件包含以下字段：
- 文件名
- 相对路径
- 是否二聚体
- 链数量
- 链ID
- 每链残基数
- 错误信息（如有）

## 注意事项

- 脚本使用BioPython库解析PDB文件，可以处理大多数标准PDB格式
- 某些特殊格式的PDB文件可能导致解析错误
- 空链ID会被自动忽略
- 脚本会统计每条链的残基数量，有助于识别有意义的链 