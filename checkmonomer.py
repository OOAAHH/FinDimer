from Bio import PDB
import os
import argparse
import csv
from tabulate import tabulate
import pandas as pd
from tqdm import tqdm
import numpy as np
from collections import defaultdict
import math

def check_dimer(pdb_file):
    """检查PDB文件是否包含多条链"""
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)
        
        # 获取所有链ID
        chains = list(structure.get_chains())
        chain_ids = [chain.id for chain in chains]
        
        # 移除空格作为链ID的情况（有些单链PDB文件使用空格作为链ID）
        chain_ids = [chain_id for chain_id in chain_ids if chain_id.strip()]
        
        # 统计每条链的残基数量，用于判断是否是有意义的链
        residue_counts = {}
        for chain in chains:
            if chain.id.strip():  # 排除空链ID
                residue_counts[chain.id] = len(list(chain.get_residues()))
        
        # 分析链特征
        chain_features = analyze_chain_features(structure)
        
        is_dimer = len(set(chain_ids)) > 1
        
        return {
            "filename": os.path.basename(pdb_file),
            "is_dimer": is_dimer,
            "chain_count": len(set(chain_ids)),
            "chain_ids": list(set(chain_ids)),
            "residue_counts": residue_counts,
            "chain_features": chain_features
        }
    except Exception as e:
        return {
            "filename": os.path.basename(pdb_file),
            "error": str(e),
            "is_dimer": False,
            "chain_count": 0,
            "chain_ids": []
        }

# 批量检测目录中所有PDB文件
def batch_check_dimers(directory, recursive=False):
    results = []
    
    # 获取所有PDB文件
    pdb_files = []
    if recursive:
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.lower().endswith(".pdb"):  # 兼容大小写扩展名
                    pdb_files.append(os.path.join(root, file))
    else:
        pdb_files = [os.path.join(directory, file) for file in os.listdir(directory) 
                    if file.lower().endswith(".pdb") and os.path.isfile(os.path.join(directory, file))]
    
    if not pdb_files:
        print(f"警告: 在目录 '{directory}' 中未找到PDB文件")
        if not recursive:
            print("提示: 如果PDB文件可能在子目录中，请使用 --recursive 参数")
        return results
    
    # 使用进度条处理文件
    for pdb_path in tqdm(pdb_files, desc="处理PDB文件"):
        result = check_dimer(pdb_path)
        # 添加相对路径信息
        result["relative_path"] = os.path.relpath(pdb_path, directory)
        results.append(result)
    
    return results

def summarize_results(results):
    """汇总检测结果"""
    total_files = len(results)
    dimer_count = sum(1 for r in results if r.get("is_dimer", False))
    monomer_count = total_files - dimer_count
    error_count = sum(1 for r in results if "error" in r)
    
    chain_counts = {}
    for result in results:
        count = result.get("chain_count", 0)
        chain_counts[count] = chain_counts.get(count, 0) + 1
    
    # 统计模块分布
    module_stats = defaultdict(int)
    for result in results:
        if "chain_features" in result and "modules" in result["chain_features"]:
            for chain_id, modules in result["chain_features"]["modules"].items():
                for module_type, count in modules.items():
                    module_stats[module_type] += count
    
    summary = {
        "总文件数": total_files,
        "二聚体数量": dimer_count,
        "单体数量": monomer_count,
        "处理出错数": error_count,
        "二聚体比例": f"{dimer_count/total_files*100:.1f}%" if total_files > 0 else "0%",
        "链数分布": chain_counts,
        "结构模块统计": dict(module_stats)
    }
    
    return summary

def export_to_csv(results, output_file):
    """导出结果到CSV文件"""
    # 准备CSV数据
    csv_data = []
    for result in results:
        # 处理复杂字段
        chain_ids_str = ",".join(result.get("chain_ids", []))
        
        # 如果有残基计数信息，格式化为字符串
        residue_info = ""
        if "residue_counts" in result and result["residue_counts"]:
            residue_info = "; ".join([f"{chain}: {count}" for chain, count in result["residue_counts"].items()])
        
        row = {
            "文件名": result["filename"],
            "相对路径": result.get("relative_path", ""),
            "是否二聚体": result.get("is_dimer", False),
            "链数量": result.get("chain_count", 0),
            "链ID": chain_ids_str,
            "每链残基数": residue_info
        }
        
        # 如果有链特征信息，添加到行
        if "chain_features" in result:
            features = result["chain_features"]
            if "nucleotide_counts" in features:
                row["核苷酸数"] = "; ".join([f"{chain}: {count}" for chain, count in features["nucleotide_counts"].items()])
            if "chain_similarity" in features:
                row["链间相似度"] = features["chain_similarity"]
            if "contact_area" in features:
                row["接触面积"] = features["contact_area"]
            
            # 添加模块信息
            if "modules" in features:
                modules_info = []
                for chain_id, modules in features["modules"].items():
                    module_strs = [f"{module_type}: {count}" for module_type, count in modules.items()]
                    modules_info.append(f"{chain_id}: {{{', '.join(module_strs)}}}")
                row["结构模块"] = "; ".join(modules_info)
            
            # 添加GC含量
            if "gc_content" in features:
                row["GC含量"] = "; ".join([f"{chain}: {content:.1f}%" for chain, content in features["gc_content"].items()])
        
        # 如果有错误信息，添加到行
        if "error" in result:
            row["错误信息"] = result["error"]
        
        csv_data.append(row)
    
    # 创建DataFrame并保存为CSV
    df = pd.DataFrame(csv_data)
    df.to_csv(output_file, index=False, encoding='utf-8-sig')  # 使用带BOM的UTF-8编码以支持中文
    
    return output_file

def calculate_distance(atom1, atom2):
    """计算两个原子之间的距离"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1.coord, atom2.coord)))

def detect_base_pairs(residues):
    """检测碱基对，返回可能的配对残基索引列表"""
    base_pairs = []
    
    # RNA碱基对的典型距离范围 (Å)
    MAX_BASE_PAIR_DISTANCE = 12.0  # 保守估计
    
    # 检查每对残基
    for i in range(len(residues)):
        for j in range(i + 3, len(residues)):  # 至少间隔3个残基才考虑可能形成碱基对
            res_i = residues[i]
            res_j = residues[j]
            
            # 检查是否是RNA残基
            if res_i.get_resname().strip() not in ['A', 'U', 'G', 'C']:
                continue
            if res_j.get_resname().strip() not in ['A', 'U', 'G', 'C']:
                continue
            
            # 尝试获取N1/N9原子（碱基中心）
            try:
                # 对于嘌呤 (A,G)，使用N9；对于嘧啶 (C,U)，使用N1
                atom_i_name = 'N9' if res_i.get_resname() in ['A', 'G'] else 'N1'
                atom_j_name = 'N9' if res_j.get_resname() in ['A', 'G'] else 'N1'
                
                atom_i = res_i[atom_i_name]
                atom_j = res_j[atom_j_name]
                
                distance = calculate_distance(atom_i, atom_j)
                
                # 如果原子间距离在合理范围内，可能形成碱基对
                if distance < MAX_BASE_PAIR_DISTANCE:
                    base_pairs.append((i, j))
            except KeyError:
                # 如果找不到关键原子，跳过
                continue
    
    return base_pairs

def identify_rna_modules(residues, base_pairs):
    """识别RNA结构模块"""
    modules = defaultdict(int)
    
    if len(residues) < 3:
        return dict(modules)
    
    # 根据碱基对构建二级结构图
    paired = [False] * len(residues)
    pairing = [-1] * len(residues) 
    
    for i, j in base_pairs:
        paired[i] = paired[j] = True
        pairing[i] = j
        pairing[j] = i
    
    # 识别茎结构 (至少连续2个碱基对)
    stems = []
    i = 0
    while i < len(residues) - 3:
        if paired[i] and pairing[i] > i:
            j = pairing[i]
            # 检查是否有连续碱基对
            stem_length = 1
            k = 1
            while (i + k < j - k and 
                   i + k < len(residues) and 
                   paired[i + k] and 
                   pairing[i + k] == j - k):
                stem_length += 1
                k += 1
            
            if stem_length >= 2:  # 至少2个碱基对才算茎
                stems.append((i, j, stem_length))
                modules["茎结构"] += 1
                i = i + stem_length
            else:
                i += 1
        else:
            i += 1
    
    # 识别发夹环 (茎结构后跟未配对区域)
    for stem_start, stem_end, stem_length in stems:
        if stem_start + stem_length < stem_end - stem_length - 2:
            loop_size = stem_end - stem_length - (stem_start + stem_length) + 1
            if 3 <= loop_size <= 8:  # 发夹环通常有3-8个核苷酸
                modules["发夹环"] += 1
    
    # 识别内部环 (两个茎之间的未配对区域)
    for i in range(len(stems) - 1):
        curr_stem_end = stems[i][0] + stems[i][2]
        next_stem_start = stems[i+1][0]
        if curr_stem_end < next_stem_start:
            gap = next_stem_start - curr_stem_end
            if 1 <= gap <= 5:  # 内部环通常较小
                modules["内部环"] += 1
    
    # 识别凸环 (茎一侧有未配对核苷酸)
    for i in range(len(pairing) - 1):
        if paired[i] and not paired[i+1] and i+2 < len(pairing) and paired[i+2]:
            modules["凸环"] += 1
    
    # 近似统计分支数（通过计算分叉点）
    branch_points = 0
    for i in range(len(paired) - 2):
        if paired[i] and paired[i+1] and paired[i+2]:
            if pairing[i] - pairing[i+1] != 1 or pairing[i+1] - pairing[i+2] != 1:
                branch_points += 1
    
    if branch_points > 0:
        modules["多分支结构"] = branch_points
    
    # 返回模块统计
    return dict(modules)

def calculate_gc_content(residues):
    """计算GC含量"""
    gc_count = 0
    total_count = 0
    
    for residue in residues:
        res_name = residue.get_resname().strip()
        if res_name in ['G', 'C', 'DG', 'DC', 'GUA', 'CYT']:
            gc_count += 1
            total_count += 1
        elif res_name in ['A', 'U', 'DA', 'DT', 'ADE', 'THY']:
            total_count += 1
    
    if total_count == 0:
        return 0.0
    
    return (gc_count / total_count) * 100.0

def analyze_chain_features(structure):
    """分析PDB结构中各链的特征"""
    
    # 统计每条链的核苷酸数
    nucleotide_counts = {}
    
    # GC含量
    gc_content = {}
    
    # 识别RNA结构模块
    modules = {}
    
    # 计算链间接触面积 (简化版)
    contact_area = 0
    
    # 链间序列相似度
    chain_similarity = 0
    
    try:
        chains = [chain for chain in structure.get_chains() if chain.id.strip()]
        
        # 统计每条链的核苷酸数和结构模块
        for chain in chains:
            nucleotide_counts[chain.id] = 0
            rna_residues = []
            
            # 只计算RNA残基 (A, U, G, C)
            for residue in chain.get_residues():
                res_name = residue.get_resname().strip()
                if res_name in ['A', 'U', 'G', 'C', 'DA', 'DT', 'DG', 'DC', 'ADE', 'THY', 'GUA', 'CYT']:
                    nucleotide_counts[chain.id] += 1
                    rna_residues.append(residue)
            
            # 计算GC含量
            gc_content[chain.id] = calculate_gc_content(rna_residues)
            
            # 如果链足够长，分析结构模块
            if len(rna_residues) >= 5:
                # 检测碱基对
                base_pairs = detect_base_pairs(rna_residues)
                
                # 识别结构模块
                modules[chain.id] = identify_rna_modules(rna_residues, base_pairs)
        
        # 简化版链间接触计算
        # 实际应用中需要更复杂的算法
        if len(chains) > 1:
            # 简单示例: 计算第一条链和第二条链的接触
            # 在实际应用中，应该为每对链计算 SASA (溶剂可及表面积) 的变化
            contact_area = 0  # 占位值
        
        # 简化版链间序列相似度计算
        if len(chains) > 1 and len(chains) == 2:
            # 计算两条链的长度比例，作为相似度的粗略估计
            lengths = [nucleotide_counts[chain.id] for chain in chains if chain.id in nucleotide_counts]
            if len(lengths) == 2 and min(lengths) > 0:
                chain_similarity = min(lengths) / max(lengths) * 100.0
            
    except Exception as e:
        print(f"分析链特征时出错: {str(e)}")
    
    return {
        "nucleotide_counts": nucleotide_counts,
        "gc_content": gc_content,
        "modules": modules,
        "contact_area": contact_area,
        "chain_similarity": chain_similarity
    }

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='检测PDB文件中是否存在多条链（二聚体）')
    parser.add_argument('directory', help='PDB文件所在目录')
    parser.add_argument('--recursive', '-r', action='store_true', help='是否递归搜索子目录')
    parser.add_argument('--output', '-o', help='导出结果到CSV文件')
    parser.add_argument('--simple', '-s', action='store_true', help='仅显示简要统计信息')
    
    args = parser.parse_args()
    
    # 检查目录是否存在
    if not os.path.isdir(args.directory):
        print(f"错误: 目录 '{args.directory}' 不存在")
        return
    
    print(f"开始检测目录: {args.directory}" + (" (递归搜索)" if args.recursive else ""))
    results = batch_check_dimers(args.directory, args.recursive)
    
    if not results:
        print("未找到PDB文件或处理过程中出现错误")
        return
    
    # 汇总结果
    summary = summarize_results(results)
    
    # 显示汇总信息
    print("\n=== 检测结果汇总 ===")
    for key, value in summary.items():
        if key not in ["链数分布", "结构模块统计"]:
            print(f"{key}: {value}")
    
    print("\n链数分布:")
    for count, num in sorted(summary["链数分布"].items()):
        print(f"  {count}条链: {num}个文件 ({num/summary['总文件数']*100:.1f}%)")
    
    if "结构模块统计" in summary and summary["结构模块统计"]:
        print("\n结构模块统计:")
        for module_type, count in summary["结构模块统计"].items():
            print(f"  {module_type}: {count}个")
    
    # 如果不是简要模式，显示详细信息
    if not args.simple:
        print("\n=== 详细信息 ===")
        # 提取需要显示的字段
        table_data = []
        for result in results:
            if "error" in result:
                table_data.append([
                    result["filename"], 
                    "处理出错", 
                    0, 
                    "", 
                    result["error"]
                ])
            else:
                # 获取每条链的核苷酸数，用于更准确判断是否为生物学意义的二聚体
                nuc_counts = []
                if "chain_features" in result and "nucleotide_counts" in result["chain_features"]:
                    nuc_counts = list(result["chain_features"]["nucleotide_counts"].values())
                
                # 如果有多条链，并且每条链的核苷酸数都不少于10，则更可能是真正的二聚体
                meaningful_dimer = result.get("is_dimer", False) and all(count >= 10 for count in nuc_counts)
                
                # 生成模块信息字符串
                module_info = ""
                if "chain_features" in result and "modules" in result["chain_features"]:
                    modules_strs = []
                    for chain_id, mods in result["chain_features"]["modules"].items():
                        if mods:
                            mods_str = ", ".join([f"{k}:{v}" for k, v in mods.items()])
                            modules_strs.append(f"{chain_id}:[{mods_str}]")
                    module_info = "; ".join(modules_strs)
                
                table_data.append([
                    result["filename"], 
                    "是 (有意义)" if meaningful_dimer else ("是" if result["is_dimer"] else "否"), 
                    result["chain_count"], 
                    ", ".join(result["chain_ids"]),
                    ", ".join([f"{c}:{n}" for c, n in result.get("residue_counts", {}).items()]) if result.get("residue_counts") else "",
                    module_info
                ])
        
        # 使用tabulate生成表格
        print(tabulate(
            table_data, 
            headers=["文件名", "是否二聚体", "链数量", "链ID", "每链残基数", "结构模块"],
            tablefmt="grid"
        ))
    
    # 如果指定了输出文件，导出结果
    if args.output:
        output_file = export_to_csv(results, args.output)
        print(f"\n结果已导出到: {output_file}")
    
    # 打印二聚体文件列表
    dimers = [r["filename"] for r in results if r.get("is_dimer", False)]
    if dimers:
        print("\n检测到的二聚体文件:")
        for dimer in dimers:
            print(f"  - {dimer}")
    
if __name__ == "__main__":
    # 直接调用main函数
    main()