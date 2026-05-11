#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orca批量反应计算工作流脚本（csv版）
从csv读取反应对并并行进行计算，反应物/产物结构来自 save_dir/result_1/R.xyz 与 P.xyz
"""

import os
import sys
import csv
import shutil
import subprocess
import multiprocessing as mp
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import json
import fcntl  # 用于文件锁
import re  # 用于清理名称字符串中的特殊字符
import numpy as np
# import networkx as nx
# from networkx.algorithms import isomorphism as iso
# from pymatgen.core import Molecule
# from pymatgen.analysis.graphs import MoleculeGraph
# from pymatgen.analysis.local_env import OpenBabelNN
from openbabel import pybel  # 延迟导入，避免环境无插件时报错

# 常数
HARTREE_TO_KCAL = 627.509

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('orca_batch_workflow.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class OrcaBatchWorkflow:
    """orca批量反应计算工作流类（仅csv输入）"""
    
    def __init__(self, csv_path: str, output_base_dir: str, 
                 orca_exec: str = "orca", max_workers: Optional[int] = None,
                 force_recalculate: bool = False, index: Optional[int] = None,
                 default_charge: int = 0, default_multiplicity: int = 1,
                 save_csv_path: Optional[str] = None):
        """
        初始化工作流
        
        Args:
            csv_path: 反应对csv文件路径
            output_base_dir: 输出基础目录
            orca_exec: orca可执行文件路径
            max_workers: 最大并行进程数
            force_recalculate: 是否强制重新计算（忽略已有结果）
            index: 只处理csv中的第index行（从0开始）
            default_charge: 默认电荷
            default_multiplicity: 默认自旋多重度
            save_csv_path: 保存结果csv路径
        """
        self.csv_path = Path(csv_path)
        # 为避免输出冲突，index存在时为输出目录追加子目录
        self.index = index
        if index is not None:
            self.output_base_dir = Path(output_base_dir) / f"idx_{index}"
        else:
            self.output_base_dir = Path(output_base_dir)
        self.orca_exec = orca_exec
        self.max_workers = max_workers or mp.cpu_count()
        self.force_recalculate = force_recalculate
        self.default_charge = default_charge
        self.default_multiplicity = default_multiplicity
        self.reaction_rows: List[Dict[str, Any]] = []
        self.reaction_meta: Dict[str, Dict[str, Any]] = {}
        self.missing_inputs: List[Dict[str, Any]] = []
        self.save_csv_path: Optional[Path] = Path(save_csv_path) if save_csv_path else None
        
        if not self.csv_path.exists():
            raise FileNotFoundError(f"反应csv不存在: {self.csv_path}")
            
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        
        self.reaction_rows = self._load_rows_from_csv(self.csv_path)
        if self.index is not None:
            if self.index < 0 or self.index >= len(self.reaction_rows):
                raise IndexError(f"index {self.index} 超出csv长度 {len(self.reaction_rows)}")
            self.reaction_rows = [self.reaction_rows[self.index]]
        logger.info(f"初始化完成: 来自csv的反应对数={len(self.reaction_rows)}")
        logger.info(f"输出目录: {self.output_base_dir}")
        logger.info(f"最大并行数: {self.max_workers}")
        logger.info(f"强制重算: {self.force_recalculate}")

    def _load_rows_from_csv(self, csv_path: Path) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames:
                raise ValueError(f"csv缺少表头: {csv_path}")
            for idx, row in enumerate(reader):
                row["_row_index"] = idx
                rows.append(row)
        return rows

    def _normalize_index_value(self, value: Any) -> str:
        if value is None:
            return ""
        text = str(value).strip()
        if not text:
            return ""
        try:
            return str(int(float(text)))
        except Exception:
            return re.sub(r"\s+", "", text)

    def _reaction_name(self, row: Dict[str, Any], idx: int) -> str:
        """统一的反应命名：优先使用 origin_idx+aug_idx，否则使用行号"""
        origin = self._normalize_index_value(row.get("origin_idx"))
        aug = self._normalize_index_value(row.get("aug_idx"))
        if origin and aug:
            base = f"origin{origin}_aug{aug}"
        elif origin:
            base = f"origin{origin}"
        elif aug:
            base = f"aug{aug}"
        else:
            base = f"row{idx}"
        return re.sub(r"[^A-Za-z0-9._-]+", "_", base)

    def _resolve_xyz_paths(self, save_dir: Optional[str]) -> Tuple[Optional[Path], Optional[Path]]:
        if not save_dir:
            return None, None
        save_dir = str(save_dir).strip()
        if not save_dir:
            return None, None
        base = Path(save_dir) / "result_1"
        return base / "R.xyz", base / "P.xyz"

    def _row_meta(self, row: Dict[str, Any], row_index: int,
                  reactant_src: Optional[Path], product_src: Optional[Path]) -> Dict[str, Any]:
        return {
            "row_index": row_index,
            "origin_idx": row.get("origin_idx"),
            "origin_rsmi": row.get("origin_rsmi"),
            "aug_idx": row.get("aug_idx"),
            "working_dir": row.get("working_dir"),
            "save_dir": row.get("save_dir"),
            "reactant_src": str(reactant_src) if reactant_src else None,
            "product_src": str(product_src) if product_src else None,
            "original_folder": row.get("save_dir") or row.get("working_dir") or ""
        }

    def _init_result(self, reaction_name: str, meta: Dict[str, Any],
                     reaction_folder: Optional[Path] = None) -> Dict[str, Any]:
        return {
            "reaction_name": reaction_name,
            "row_index": meta.get("row_index"),
            "origin_idx": meta.get("origin_idx"),
            "origin_rsmi": meta.get("origin_rsmi"),
            "aug_idx": meta.get("aug_idx"),
            "working_dir": meta.get("working_dir"),
            "save_dir": meta.get("save_dir"),
            "reactant_src": meta.get("reactant_src"),
            "product_src": meta.get("product_src"),
            "original_folder": meta.get("original_folder"),
            "reaction_folder": reaction_folder,
            "work_dir": None,
            "success": False,
            "skipped": False,
            "reactant_energy": None,
            "product_energy": None,
            "ts_energy": None,
            "barrier_height": None,
            "reverse_barrier_height": None,
            "reaction_energy": None,
            "frequencies": [],
            "has_imaginary_freq": False,
            "imaginary_freq_count": 0,
            "is_valid_ts": False,
            "irc_analysis": {
                "performed": False,
                "success": False,
                "reactant_isomorphic": False,
                "product_isomorphic": False,
                "both_isomorphic": False
            },
            "error_message": None
        }

    def _missing_input_result(self, reaction_name: str, meta: Dict[str, Any],
                              work_dir: Path, message: str) -> Dict[str, Any]:
        result = self._init_result(reaction_name, meta)
        result["work_dir"] = work_dir
        result["success"] = False
        result["skipped"] = False
        result["error_message"] = message
        return result
    
    def is_reaction_completed(self, reaction_name: str) -> bool:
        """
        检查反应是否已经计算完成
        
        Args:
            reaction_name: 反应名称
            
        Returns:
            是否已完成
        """
        work_dir = self.output_base_dir / reaction_name
        
        # 检查必要的输出文件是否存在
        required_files = [
            work_dir / "reaction_data.json",
            work_dir / "summary.txt",
            work_dir / "reactant_opt.xyz",
            work_dir / "product_opt.xyz",
            work_dir / "ts_opt.xyz"
        ]
        
        # 所有必要文件都存在才认为已完成
        if all(f.exists() for f in required_files):
            # 额外检查JSON文件的内容是否有效
            try:
                json_file = work_dir / "reaction_data.json"
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    # 检查是否计算成功
                    if data.get('calculation_info', {}).get('successful', False):
                        return True
            except Exception as e:
                logger.warning(f"检查 {reaction_name} 的JSON文件时出错: {e}")
                return False
        
        return False
    
    def acquire_lock(self, reaction_name: str) -> Optional[object]:
        """
        为反应获取锁，防止多进程重复计算
        
        Args:
            reaction_name: 反应名称
            
        Returns:
            锁文件对象，如果无法获取则返回None
        """
        work_dir = self.output_base_dir / reaction_name
        work_dir.mkdir(parents=True, exist_ok=True)
        
        lock_file = work_dir / ".calculation.lock"
        
        try:
            # 打开锁文件
            lock_fd = open(lock_file, 'w')
            # 尝试获取非阻塞锁
            fcntl.flock(lock_fd.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            # 写入当前进程信息
            lock_fd.write(f"PID: {os.getpid()}\n")
            lock_fd.write(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            lock_fd.flush()
            return lock_fd
        except IOError:
            # 无法获取锁，说明其他进程正在计算
            logger.info(f"反应 {reaction_name} 正在被其他进程计算，跳过")
            return None
        except Exception as e:
            logger.error(f"获取锁失败 {reaction_name}: {e}")
            return None
    
    def release_lock(self, lock_fd):
        """释放锁"""
        if lock_fd:
            try:
                fcntl.flock(lock_fd.fileno(), fcntl.LOCK_UN)
                lock_fd.close()
            except Exception as e:
                logger.error(f"释放锁失败: {e}")
    
    def get_reaction_folders(self) -> List[Path]:
        """基于csv生成反应目录及xyz文件"""
        reaction_folders = []
        for idx, row in enumerate(self.reaction_rows):
            row_index = row.get("_row_index", idx)
            try:
                row_index = int(row_index)
            except Exception:
                row_index = idx

            reaction_name = self._reaction_name(row, row_index)
            reactant_src, product_src = self._resolve_xyz_paths(row.get("save_dir"))
            row_meta = self._row_meta(row, row_index, reactant_src, product_src)

            if reactant_src is None or product_src is None:
                message = "缺少save_dir，无法定位R/P.xyz"
                logger.warning(f"{reaction_name} {message}")
                missing = self._missing_input_result(reaction_name, row_meta,
                                                     self.output_base_dir / reaction_name, message)
                self.missing_inputs.append(missing)
                continue
            if not reactant_src.exists() or not product_src.exists():
                message = f"输入结构不存在: {reactant_src} / {product_src}"
                logger.warning(f"{reaction_name} {message}")
                missing = self._missing_input_result(reaction_name, row_meta,
                                                     self.output_base_dir / reaction_name, message)
                self.missing_inputs.append(missing)
                continue

            work_dir = self.output_base_dir / reaction_name
            work_dir.mkdir(parents=True, exist_ok=True)
            reactant_file = work_dir / "reactant.xyz"
            product_file = work_dir / "product.xyz"
            try:
                shutil.copy(reactant_src, reactant_file)
                shutil.copy(product_src, product_file)
            except Exception as e:
                message = f"复制输入结构失败: {e}"
                logger.warning(f"{reaction_name} {message}")
                missing = self._missing_input_result(reaction_name, row_meta, work_dir, message)
                self.missing_inputs.append(missing)
                continue

            self.reaction_meta[reaction_name] = {
                **row_meta,
                "reactant_charge": self.default_charge,
                "reactant_mult": self.default_multiplicity,
                "product_charge": self.default_charge,
                "product_mult": self.default_multiplicity,
            }
            reaction_folders.append(work_dir)
        logger.info(f"csv模式生成 {len(reaction_folders)} 个反应目录")
        return reaction_folders
    
    def create_orca_input(self, xyz_file: Path, calc_type: str, 
                         additional_keywords: str = "", reactant_file: Path = None,
                         charge: int = 0, multiplicity: int = 1) -> str:
        """
        创建orca输入文件
        
        Args:
            xyz_file: xyz文件路径
            calc_type: 计算类型 (opt, neb, ts_opt, freq)
            additional_keywords: 额外关键词
            reactant_file: 反应物文件路径（NEB计算需要）
            
        Returns:
            orca输入文件内容
        """
        # 读取xyz文件内容并解析
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        # 跳过第一行（原子数）和第二行
        coord_lines = []
        for line in lines[2:]:
            coord_lines.append(line)

        xyz_content = '\n'.join(coord_lines)
        
        # 根据计算类型设置关键词
        if calc_type == "opt":
            # keywords = "! wB97X def2-svp def2/J OPT RIJCOSX"
            keywords = "! XTB LOOSEOPT"
        elif calc_type == "neb":
            keywords = "! XTB NEB-CI"
            # keywords = "! wB97X def2-svp def2/J RIJCOSX NEB-CI"
        elif calc_type == "ts_opt":
            # keywords = "! wB97X def2-svp def2/J OptTS RIJCOSX"
            keywords = "! XTB OptTS"
        elif calc_type == "freq":
            # keywords = "! wB97X def2-svp def2/J FREQ RIJCOSX"
            keywords = "! XTB FREQ"
        else:
            raise ValueError(f"未知的计算类型: {calc_type}")
        
        # 添加额外关键词
        if additional_keywords:
            keywords += f" {additional_keywords}"
        
        # 构建输入文件
        if calc_type == "neb" and reactant_file:
            # NEB计算需要特殊的输入格式（保持1核）
            # Tol_MaxFP_I      5.e-3
            # Tol_RMSFP_I      2.5.e-3
            # Tol_MaxF_CI      5.e-4
            # Tol_RMSF_CI      2.5.e-4
            input_content = f"""%pal
nprocs 1
end

{keywords}

%neb
neb_end_xyzfile "product_opt.xyz"
NImages 7
MaxIter 200
end

%geom
MaxIter 200
end

* xyzfile {charge} {multiplicity} {reactant_file.name}
"""
        else:
            # 普通计算格式（使用8核）
            if calc_type in ["opt", "ts_opt"]:
                input_content = f"""%pal
nprocs 1
end

{keywords}

%geom
MaxIter 200
end

* xyz {charge} {multiplicity}
{xyz_content}
*

"""
            else:
                input_content = f"""%pal
nprocs 1
end

{keywords}

* xyz {charge} {multiplicity}
{xyz_content}
*

"""
        return input_content
    
    def run_orca_calculation(self, input_content: str, output_prefix: str, 
                         work_dir: Path) -> Tuple[bool, str]:
        """运行orca计算"""
        input_file = work_dir / f"{output_prefix}.inp"
        output_file = work_dir / f"{output_prefix}.runner.log"

        try:
            # 写入输入文件
            with open(input_file, 'w', encoding="utf-8") as f:
                f.write(input_content)

            # 保存当前目录，并切换到工作目录
            old_cwd = os.getcwd()
            os.chdir(work_dir)

            try:
                # 运行 ORCA
                cmd = [self.orca_exec, str(input_file.name)]
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=360000  # 100小时超时
                )
            finally:
                # 回到原始目录
                os.chdir(old_cwd)

            # 保存捕获的stdout/stderr信息
            with open(output_file, 'w', encoding="utf-8") as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
                f.write(f"\n\nReturn code: {result.returncode}\n")

            success = (result.returncode == 0)
            return success, result.stdout + result.stderr

        except subprocess.TimeoutExpired:
            logger.error(f"orca计算超时: {output_prefix}")
            return False, "计算超时"
        except Exception as e:
            logger.error(f"orca计算失败: {output_prefix}, 错误: {e}")
            return False, str(e)
    
    def extract_energy(self, output_content: str) -> Optional[float]:
        """从orca输出中提取能量"""
        lines = output_content.split('\n')
        for line in lines:
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    energy = float(line.split()[-1])
                    return energy
                except (ValueError, IndexError):
                    continue
        return None
    
    def extract_frequencies(self, output_content: str) -> List[float]:
        """从频率计算输出中提取频率"""
        frequencies = []
        lines = output_content.split('\n')
        in_freq_section = False
        
        for line in lines:
            if "VIBRATIONAL FREQUENCIES" in line:
                in_freq_section = True
                continue
            elif in_freq_section:
                if "cm**-1" in line and ":" in line:
                    try:
                        parts = line.split(":")
                        if len(parts) >= 2:
                            freq_part = parts[1].strip()
                            freq_str = freq_part.split()[0]
                            freq = float(freq_str)
                            frequencies.append(freq)
                    except (ValueError, IndexError):
                        continue
                elif line.strip() and not line.startswith(' ') and not line.startswith('\t') and not line.startswith('-'):
                    if not any(keyword in line.lower() for keyword in ['scaling', 'imaginary', 'mode', 'frequencies', 'applied']):
                        break
        
        return frequencies
    
    def count_imaginary_frequencies(self, frequencies: List[float]) -> int:
        """计算虚频数量"""
        return sum(1 for freq in frequencies if freq < 0)
    
    def is_valid_transition_state(self, frequencies: List[float]) -> bool:
        """检查是否为有效的过渡态（恰好1个虚频）"""
        imaginary_count = self.count_imaginary_frequencies(frequencies)
        return imaginary_count == 1
    
    def create_irc_input(self, ts_xyz_file: Path, work_dir: Path, charge: int, multiplicity: int) -> str:
        """创建IRC计算输入文件"""
        # 读取过渡态xyz文件
        with open(ts_xyz_file, 'r') as f:
            lines = f.readlines()
        
        # 跳过第一行（原子数）和第二行
        coord_lines = []
        for line in lines[2:]:
            coord_lines.append(line)
        
        xyz_content = '\n'.join(coord_lines)
        
        # 创建IRC输入文件（使用ωB97X/6-31G*方法，8核）
        # wB97X def2-svp def2/J RIJCOSX
        input_content = f"""%pal
   nprocs 1
end

! xTB IRC

%geom
   MaxIter 200
   Recalc_Hess 5
end

%irc
   InitHess   read
   Hess_Filename "freq.hess"
   Direction both          # 正反两个方向都跟踪IRC
   MaxIter 200             # 每个方向最多走200步
   Scale_Init_Displ 0.1    # 初始步长（步长控制参数）
   Scale_Displ_SD 0.15     # 每一步的标准位移缩放
end

* xyz {charge} {multiplicity}
{xyz_content}
*

"""
        return input_content
    
    def run_irc_calculation(self, ts_xyz_file: Path, work_dir: Path, charge: int, multiplicity: int) -> Tuple[bool, str]:
        """运行IRC计算"""
        logger.info(f"开始IRC计算")
        
        # 创建IRC输入文件
        irc_input = self.create_irc_input(ts_xyz_file, work_dir, charge, multiplicity)
        
        # 运行IRC计算
        success, output = self.run_orca_calculation(
            irc_input, "irc", work_dir
        )
        
        return success, output
    
    def extract_irc_structures(self, irc_output: str, work_dir: Path) -> Tuple[bool, Optional[Path], Optional[Path]]:
        """从IRC输出中提取反应物和产物端结构"""
        irc_reactant_path = work_dir / "irc_IRC_B.xyz"
        irc_product_path = work_dir / "irc_IRC_F.xyz"

        if irc_reactant_path.exists() and irc_product_path.exists():
            return True, irc_reactant_path, irc_product_path
        else:
            return False, None, None
        
            
    
    def _extract_structures_from_output(self, output_content: str, work_dir: Path) -> Tuple[bool, Optional[Path], Optional[Path]]:
        """从输出文本中提取IRC结构"""
        try:
            lines = output_content.split('\n')
            structures = []
            current_structure = []
            in_structure = False
            atom_count = 0
            
            for line in lines:
                line = line.strip()
                
                # 查找CARTESIAN COORDINATES部分
                if 'CARTESIAN COORDINATES' in line:
                    in_structure = True
                    current_structure = []
                    continue
                
                if in_structure:
                    if line and not line.startswith('--') and not line.startswith('='):
                        parts = line.split()
                        if len(parts) >= 4 and parts[0].isalpha():
                            try:
                                # 验证是否为坐标行
                                float(parts[1])
                                float(parts[2]) 
                                float(parts[3])
                                current_structure.append(line)
                                continue
                            except (ValueError, IndexError):
                                pass
                    
                    # 如果遇到空行或分隔符，保存当前结构
                    if (not line or line.startswith('--') or line.startswith('=')) and current_structure:
                        if len(current_structure) > 0:
                            structures.append(current_structure)
                            current_structure = []
                            in_structure = False
            
            # 保存最后一个结构
            if current_structure:
                structures.append(current_structure)
            
            if len(structures) < 2:
                logger.warning("从输出中提取的结构数量不足")
                return False, None, None
            
            # 保存第一个和最后一个结构
            reactant_path = work_dir / "irc_reactant.xyz"
            product_path = work_dir / "irc_product.xyz"
            
            # 写入反应物端结构
            with open(reactant_path, 'w') as f:
                f.write(f"{len(structures[0])}\n")
                f.write("IRC reactant structure\n")
                for line in structures[0]:
                    f.write(f"{line}\n")
            
            # 写入产物端结构
            with open(product_path, 'w') as f:
                f.write(f"{len(structures[-1])}\n")
                f.write("IRC product structure\n")
                for line in structures[-1]:
                    f.write(f"{line}\n")
            
            logger.info(f"从输出中提取IRC结构: {reactant_path.name}, {product_path.name}")
            return True, reactant_path, product_path
            
        except Exception as e:
            logger.error(f"从输出中提取结构失败: {e}")
            return False, None, None

    def canonical_smiles_from_xyz(self, xyz_file: Path) -> str:
        """尽量将 xyz 转为 canonical SMILES（含立体）；若依赖缺失则返回空"""
        try:
            mol = next(pybel.readfile("xyz", str(xyz_file)))
            return mol.write("can").strip().split()[0]
        except Exception as e:
            logger.warning(f"canonical_smiles_from_xyz 失败 {xyz_file}: {e}")
            return ""
    
    def cross_isomorphism(self, true_r: Path, true_p: Path, irc_r: Path, irc_p: Path) -> Dict[str, Any]:
        """交叉比对同构性，使用canonical SMILES（保留立体信息）"""
        smi_true_r = self.canonical_smiles_from_xyz(true_r)
        smi_true_p = self.canonical_smiles_from_xyz(true_p)
        smi_irc_r = self.canonical_smiles_from_xyz(irc_r)
        smi_irc_p = self.canonical_smiles_from_xyz(irc_p)

        if not all([smi_true_r, smi_true_p, smi_irc_r, smi_irc_p]):
            return {"performed": True, "success": False, "endpoint_match": "Error"}

        matches = (
            smi_true_r == smi_irc_r,
            smi_true_p == smi_irc_p,
            smi_true_r == smi_irc_p,
            smi_true_p == smi_irc_r,
        )
        rxn_status = "Conformational change" if smi_irc_r == smi_irc_p else "Chemical reaction"

        if matches in [(True, True, False, False), (False, False, True, True)]:
            endpoint_match = "2-end match"
        elif any(matches):
            endpoint_match = "1-end match"
        else:
            endpoint_match = "No match"
        return {
            "performed": True,
            "success": True,
            "endpoint_match": endpoint_match,
            "matches": matches,
            "rxn_status": rxn_status,
        }


    def write_xyz_from_output(self, output_content: str, xyz_path: Path, title: str) -> bool:
        """从orca输出文本中提取最后一次的笛卡尔坐标并写入xyz文件"""
        lines = output_content.split('\n')
        blocks = []
        in_block = False
        current_block = []

        for line in lines:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                in_block = True
                current_block = []
                continue
            if in_block:
                if line.strip() == '':
                    continue
                parts = line.split()
                if len(parts) >= 4 and parts[0].isalpha():
                    try:
                        _ = float(parts[1]); _ = float(parts[2]); _ = float(parts[3])
                        current_block.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                        continue
                    except Exception:
                        if current_block:
                            blocks.append(current_block)
                        in_block = False
                else:
                    if current_block:
                        blocks.append(current_block)
                        in_block = False

        if in_block and current_block:
            blocks.append(current_block)

        if not blocks:
            return False

        final_block = blocks[-1]
        try:
            with open(xyz_path, 'w') as f:
                f.write(f"{len(final_block)}\n")
                f.write(f"{title}\n")
                for sym, x, y, z in final_block:
                    f.write(f"{sym} {x:.6f} {y:.6f} {z:.6f}\n")
            logger.info(f"文件已保存到{xyz_path}")
            return True
        except Exception:
            return False
    
    def write_xyz_from_G_array(self, elements: List[str], coords: Any, xyz_path: Path, title: str = ""):
        """根据元素列表和坐标数组写xyz文件"""
        if elements is None or coords is None or len(elements) == 0:
            raise ValueError("缺少元素或坐标")
        coords_arr = np.asarray(coords, dtype=float)
        if coords_arr.ndim != 2 or coords_arr.shape[0] == 0:
            raise ValueError("坐标为空或格式错误")
        if coords_arr.shape[0] != len(elements):
            raise ValueError("坐标与元素数量不匹配")
        valid_mask = np.isfinite(coords_arr[:, :3]).all(axis=1)
        if not np.all(valid_mask):
            coords_arr = coords_arr[valid_mask]
            elements = [elem for elem, keep in zip(elements, valid_mask) if keep]
        if coords_arr.shape[0] == 0 or len(elements) == 0:
            raise ValueError("所有坐标均为无效值")
        with open(xyz_path, "w") as f:
            f.write(f"{len(elements)}\n")
            f.write(f"{title}\n")
            for elem, (x, y, z) in zip(elements, coords_arr[:, :3]):
                f.write(f"{elem} {x:.6f} {y:.6f} {z:.6f}\n")
        return xyz_path
    
    def xyz_to_array(self, xyz_path: Path) -> Optional[np.ndarray]:
        """读取xyz并返回坐标数组"""
        if not xyz_path.exists():
            return None
        try:
            with open(xyz_path, "r") as f:
                lines = f.readlines()
            n_atoms = int(lines[0].strip())
            coords = []
            for line in lines[2:2 + n_atoms]:
                parts = line.strip().split()
                if len(parts) >= 4:
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            if len(coords) != n_atoms:
                return None
            return np.array(coords, dtype=float)
        except Exception:
            return None
    
    def read_xyz_file(self, xyz_path: Path) -> Optional[str]:
        """读取xyz文件内容"""
        try:
            if xyz_path.exists():
                with open(xyz_path, 'r', encoding='utf-8') as f:
                    return f.read()
            return None
        except Exception as e:
            logger.warning(f"读取xyz文件失败 {xyz_path}: {e}")
            return None
    
    def process_single_reaction(self, reaction_folder: Path) -> Dict:
        """
        处理单个反应文件夹的完整工作流（带锁和状态检查）
        
        Args:
            reaction_folder: 反应文件夹路径（包含reactant.xyz和product.xyz）
            
        Returns:
            处理结果字典
        """
        # 读取反应物和产物文件
        reactant_file = reaction_folder / "reactant.xyz"
        product_file = reaction_folder / "product.xyz"
        
        reaction_name = reaction_folder.name
        meta = self.reaction_meta.get(reaction_name, {})
        r_charge = meta.get("reactant_charge", 0)
        r_mult = meta.get("reactant_mult", 1)
        p_charge = meta.get("product_charge", 0)
        p_mult = meta.get("product_mult", 1)
        base_result = self._init_result(reaction_name, meta, reaction_folder)

        if not reactant_file.exists() or not product_file.exists():
            base_result["error_message"] = "反应物或产物结构文件不存在"
            return base_result
        
        # 检查是否已经计算完成
        if not self.force_recalculate and self.is_reaction_completed(reaction_name):
            logger.info(f"反应 {reaction_name} 已完成，跳过计算")
            # 读取已有结果
            work_dir = self.output_base_dir / reaction_name
            json_file = work_dir / "reaction_data.json"
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                barrier_analysis = data.get('barrier_analysis', {})
                energies = data.get('energies', {})
                freq_analysis = data.get('frequency_analysis', {})
                irc_analysis = data.get('irc_analysis', base_result.get("irc_analysis"))
                base_result.update({
                    'success': True,
                    'skipped': True,
                    'work_dir': work_dir,
                    'reactant_energy': energies.get('reactant_energy_hartree'),
                    'product_energy': energies.get('product_energy_hartree'),
                    'ts_energy': energies.get('ts_energy_hartree'),
                    'barrier_height': barrier_analysis.get('barrier_height_hartree'),
                    'reverse_barrier_height': barrier_analysis.get('reverse_barrier_height_hartree'),
                    'reaction_energy': barrier_analysis.get('reaction_energy_hartree'),
                    'frequencies': freq_analysis.get('all_frequencies', []),
                    'has_imaginary_freq': freq_analysis.get('has_imaginary_freq', False),
                    'imaginary_freq_count': freq_analysis.get('imaginary_frequencies_count', 0),
                    'is_valid_ts': freq_analysis.get('is_valid_transition_state', False),
                    'irc_analysis': irc_analysis
                })
                return base_result
            except Exception as e:
                logger.warning(f"读取已有结果失败: {e}")
        
        # 尝试获取锁
        lock_fd = self.acquire_lock(reaction_name)
        if lock_fd is None:
            # 无法获取锁，跳过
            base_result.update({
                'success': False,
                'skipped': True,
                'error_message': '正在被其他进程计算'
            })
            return base_result
        
        try:
            logger.info(f"开始处理反应: {reaction_name} (原始文件夹: {reaction_folder.name})")
            
            # 创建输出工作目录
            work_dir = self.output_base_dir / reaction_name
            work_dir.mkdir(exist_ok=True)
            
            result = base_result
            result.update({
                'reactant_file': reactant_file,
                'product_file': product_file,
                'ts_file': None,
                'work_dir': work_dir,
                'success': False,
                'skipped': False
            })
            
            # 执行计算流程（与原代码相同）
            # 步骤1: 反应物优化
            logger.info(f"步骤1: 优化反应物 {reaction_name}")
            reactant_input = self.create_orca_input(reactant_file, "opt", charge=r_charge, multiplicity=r_mult)
            success, output = self.run_orca_calculation(
                reactant_input, "reactant_opt", work_dir
            )
            if not success:
                result['error_message'] = f"反应物优化失败"
                return result
            
            result['reactant_energy'] = self.extract_energy(output)
            self.write_xyz_from_output(output, work_dir / "reactant_opt.xyz", f"reactant optimized: {reaction_name}")
            
            # 步骤2: 产物优化
            logger.info(f"步骤2: 优化产物 {reaction_name}")
            product_input = self.create_orca_input(product_file, "opt", charge=p_charge, multiplicity=p_mult)
            success, output = self.run_orca_calculation(
                product_input, "product_opt", work_dir
            )
            if not success:
                result['error_message'] = f"产物优化失败"
                return result
            
            result['product_energy'] = self.extract_energy(output)
            self.write_xyz_from_output(output, work_dir / "product_opt.xyz", f"product optimized: {reaction_name}")
            
            # 步骤3: CI-NEB计算
            logger.info(f"步骤3: CI-NEB计算 {reaction_name}")
            reactant_neb_file = work_dir / "reactant_opt.xyz"
            product_neb_file = work_dir / "product_opt.xyz"
            
            neb_input = self.create_orca_input(product_neb_file, "neb", reactant_file=reactant_neb_file, charge=r_charge, multiplicity=r_mult)
            success, output = self.run_orca_calculation(
                neb_input, "neb", work_dir
            )
            if not success:
                result['error_message'] = f"CI-NEB计算失败"
                return result
            
            # NEB 收敛后直接使用收敛结构作为 TS 初猜
            converged_xyz = work_dir / "neb_NEB-CI_converged.xyz"
            ts_file = work_dir / "ts_ini.xyz"
            if converged_xyz.exists():
                shutil.copy(converged_xyz, ts_file)
            else:
                logger.warning("CI-NEB 收敛失败")
            
            # 步骤4: 过渡态优化
            logger.info(f"步骤4: 过渡态优化 {reaction_name}")
            if not ts_file.exists():
                result['error_message'] = f"未找到TS初始结构"
                return result
                
            ts_input = self.create_orca_input(ts_file, "ts_opt", charge=r_charge, multiplicity=r_mult)
            success, output = self.run_orca_calculation(
                ts_input, "ts_opt", work_dir
            )

            if not success:
                result['error_message'] = f"过渡态优化失败"
                return result
            
            result['ts_energy'] = self.extract_energy(output)
            self.write_xyz_from_output(output, work_dir / "ts_opt.xyz", f"TS optimized: {reaction_name}")
            
            # 步骤5: 频率计算
            logger.info(f"步骤5: 频率计算 {reaction_name}")
            ts_opt_file = work_dir / "ts_opt.xyz"
            result['ts_file']=ts_opt_file
            freq_input = self.create_orca_input(ts_opt_file, "freq", charge=r_charge, multiplicity=r_mult)
            success, output = self.run_orca_calculation(
                freq_input, "freq", work_dir
            )
            if not success:
                result['error_message'] = f"频率计算失败"
                return result
            
            # 提取频率
            result['frequencies'] = self.extract_frequencies(output)
            result['has_imaginary_freq'] = any(f < 0 for f in result['frequencies'])
            result['imaginary_freq_count'] = self.count_imaginary_frequencies(result['frequencies'])
            result['is_valid_ts'] = self.is_valid_transition_state(result['frequencies'])
            
            # 步骤6: 过渡态筛选和IRC计算
            result['irc_success'] = False
            result['irc_reactant_path'] = None
            result['irc_product_path'] = None
            result['reactant_isomorphic'] = False
            result['product_isomorphic'] = False
            result['irc_analysis'] = {
                'performed': False,
                'success': False,
                'reactant_isomorphic': False,
                'product_isomorphic': False,
                'both_isomorphic': False
            }
            
            if result['is_valid_ts']:
                logger.info(f"步骤6: 过渡态有效（虚频数=1），开始IRC计算 {reaction_name}")
                
                # 运行IRC计算
                irc_success, irc_output = self.run_irc_calculation(ts_opt_file, work_dir, r_charge, r_mult)
                result['irc_success'] = irc_success
                
                if irc_success:
                    # 提取IRC结构
                    irc_extract_success, irc_reactant_path, irc_product_path = self.extract_irc_structures(
                        irc_output, work_dir
                    )
                    
                    if irc_extract_success and irc_reactant_path and irc_product_path:
                        result['irc_reactant_path'] = irc_reactant_path
                        result['irc_product_path'] = irc_product_path
                        # 先对IRC端点再优化
                        logger.info("步骤6.1: IRC端点再优化用于同构比对")
                        irc_r_opt = work_dir / "irc_reactant_opt.xyz"
                        irc_p_opt = work_dir / "irc_product_opt.xyz"
                        irc_r_input = self.create_orca_input(irc_reactant_path, "opt", charge=r_charge, multiplicity=r_mult)
                        irc_p_input = self.create_orca_input(irc_product_path, "opt", charge=p_charge, multiplicity=p_mult)
                        irc_r_success, irc_r_out = self.run_orca_calculation(irc_r_input, "irc_reactant_opt", work_dir)
                        if irc_r_success:
                            self.write_xyz_from_output(irc_r_out, irc_r_opt, "irc reactant reopt")
                        irc_p_success, irc_p_out = self.run_orca_calculation(irc_p_input, "irc_product_opt", work_dir)
                        if irc_p_success:
                            self.write_xyz_from_output(irc_p_out, irc_p_opt, "irc product reopt")
                        result["irc_r_energy"] = self.extract_energy(irc_r_out) if irc_r_success else None
                        result["irc_p_energy"] = self.extract_energy(irc_p_out) if irc_p_success else None

                        # 步骤7: 同构性检查（交叉比对）
                        logger.info(f"步骤7: 同构性检查 {reaction_name}")
                        if irc_r_success and irc_p_success:
                            cross_res = self.cross_isomorphism(
                                work_dir / "reactant_opt.xyz",
                                work_dir / "product_opt.xyz",
                                irc_r_opt,
                                irc_p_opt,
                            )
                            result['irc_analysis'] = cross_res
                            logger.info(f"同构性检查完成，匹配: {cross_res.get('endpoint_match')}")
                        else:
                            result['irc_analysis']['performed'] = True
                            result['irc_analysis']['success'] = False
                    else:
                        logger.warning(f"IRC结构提取失败 {reaction_name}")
                        result['irc_analysis']['performed'] = True
                        result['irc_analysis']['success'] = False
                else:
                    logger.warning(f"IRC计算失败 {reaction_name}")
                    result['irc_analysis']['performed'] = True
                    result['irc_analysis']['success'] = False
            else:
                logger.warning(f"过渡态无效（虚频数={result['imaginary_freq_count']}），跳过IRC计算 {reaction_name}")
                result['irc_analysis']['performed'] = False
            
            # 计算能垒和反应能
            if result['reactant_energy'] and result['ts_energy']:
                result['barrier_height'] = result['ts_energy'] - result['reactant_energy']
            if result['product_energy'] and result['ts_energy']:
                result['reverse_barrier_height'] = result['ts_energy'] - result['product_energy']
            if result['reactant_energy'] and result['product_energy']:
                result['reaction_energy'] = result['product_energy'] - result['reactant_energy']
            
            # 保存结果文件
            result['success'] = True
            self.save_result_files(result)

            logger.info(f"反应 {reaction_name} 处理完成")
            
        except Exception as e:
            result['error_message'] = f"处理过程中发生错误: {str(e)}"
            logger.error(f"处理反应 {reaction_name} 时发生错误: {e}")
        finally:
            # 释放锁
            self.release_lock(lock_fd)
        
        return result
    
    def save_result_files(self, result: Dict):
        """保存结果文件"""
        work_dir = result['work_dir']
        
        # 保存能量和频率信息
        summary_file = work_dir / "summary.txt"
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(f"反应: {result['reaction_name']}\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("能量信息 (Hartree):\n")
            f.write("-" * 30 + "\n")
            f.write(f"反应物能量: {result['reactant_energy']:.6f}\n")
            f.write(f"产物能量: {result['product_energy']:.6f}\n")
            f.write(f"过渡态能量: {result['ts_energy']:.6f}\n\n")
            
            f.write("能量信息 (kcal/mol):\n")
            f.write("-" * 30 + "\n")
            f.write(f"反应物能量: {result['reactant_energy'] * HARTREE_TO_KCAL:.2f}\n")
            f.write(f"产物能量: {result['product_energy'] * HARTREE_TO_KCAL:.2f}\n")
            f.write(f"过渡态能量: {result['ts_energy'] * HARTREE_TO_KCAL:.2f}\n\n")
            
            f.write("反应分析:\n")
            f.write("-" * 30 + "\n")
            f.write(f"反应能垒: {result['barrier_height']:.6f} Hartree ({result['barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
            f.write(f"逆反应能垒: {result['reverse_barrier_height']:.6f} Hartree ({result['reverse_barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
            f.write(f"反应能: {result['reaction_energy']:.6f} Hartree ({result['reaction_energy'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n\n")
            
            f.write("频率分析:\n")
            f.write("-" * 30 + "\n")
            imaginary_count = sum(1 for freq in result['frequencies'] if freq < 0)
            f.write(f"虚频数量: {imaginary_count}\n")
            f.write(f"有虚频: {'是' if result['has_imaginary_freq'] else '否'}\n")
            f.write(f"有效过渡态: {'是' if result.get('is_valid_ts', False) else '否'}\n")
            if result['frequencies']:
                f.write(f"所有频率: {', '.join([f'{freq:.2f}' for freq in result['frequencies'][:10]])}...\n")
            
            f.write("\nIRC分析:\n")
            f.write("-" * 30 + "\n")
            irc_analysis = result.get('irc_analysis', {})
            f.write(f"IRC计算执行: {'是' if irc_analysis.get('performed', False) else '否'}\n")
            if irc_analysis.get('performed', False):
                f.write(f"IRC计算成功: {'是' if irc_analysis.get('success', False) else '否'}\n")
                f.write(f"反应物同构: {'是' if irc_analysis.get('reactant_isomorphic', False) else '否'}\n")
                f.write(f"产物同构: {'是' if irc_analysis.get('product_isomorphic', False) else '否'}\n")
                f.write(f"两端都同构: {'是' if irc_analysis.get('both_isomorphic', False) else '否'}\n")
        
        # 读取xyz文件内容
        reactant_xyz_content = self.read_xyz_file(work_dir / "reactant_opt.xyz")
        product_xyz_content = self.read_xyz_file(work_dir / "product_opt.xyz")
        ts_xyz_content = self.read_xyz_file(work_dir / "ts_opt.xyz")
        
        # 读取IRC结构内容（如果存在）
        irc_reactant_xyz_content = None
        irc_product_xyz_content = None
        if result.get('irc_reactant_path') and result['irc_reactant_path'].exists():
            irc_reactant_xyz_content = self.read_xyz_file(result['irc_reactant_path'])
        if result.get('irc_product_path') and result['irc_product_path'].exists():
            irc_product_xyz_content = self.read_xyz_file(result['irc_product_path'])
        
        # 计算虚频个数
        imaginary_freq_count = sum(1 for freq in result['frequencies'] if freq < 0)
        
        # 保存JSON格式数据（包含xyz内容和IRC分析）
        json_data = {
            "reaction_name": result['reaction_name'],
            "original_folder": result.get('original_folder', ''),
            "energies": {
                "reactant_energy_hartree": result['reactant_energy'],
                "product_energy_hartree": result['product_energy'],
                "ts_energy_hartree": result['ts_energy'],
                "reactant_energy_kcal_mol": result['reactant_energy'] * HARTREE_TO_KCAL if result['reactant_energy'] else None,
                "product_energy_kcal_mol": result['product_energy'] * HARTREE_TO_KCAL if result['product_energy'] else None,
                "ts_energy_kcal_mol": result['ts_energy'] * HARTREE_TO_KCAL if result['ts_energy'] else None
            },
            "barrier_analysis": {
                "barrier_height_hartree": result['barrier_height'],
                "barrier_height_kcal_mol": result['barrier_height'] * HARTREE_TO_KCAL if result['barrier_height'] else None,
                "reverse_barrier_height_hartree": result['reverse_barrier_height'],
                "reverse_barrier_height_kcal_mol": result['reverse_barrier_height'] * HARTREE_TO_KCAL if result['reverse_barrier_height'] else None,
                "reaction_energy_hartree": result['reaction_energy'],
                "reaction_energy_kcal_mol": result['reaction_energy'] * HARTREE_TO_KCAL if result['reaction_energy'] else None
            },
            "frequency_analysis": {
                "imaginary_frequencies_count": imaginary_freq_count,
                "has_imaginary_freq": result['has_imaginary_freq'],
                "is_valid_transition_state": result.get('is_valid_ts', False),
                "all_frequencies": result['frequencies']
            },
            "irc_analysis": result.get('irc_analysis', {
                'performed': False,
                'success': False,
                'reactant_isomorphic': False,
                'product_isomorphic': False,
                'both_isomorphic': False
            }),
            "xyz_structures": {
                "reactant_xyz": reactant_xyz_content,
                "product_xyz": product_xyz_content,
                "ts_xyz": ts_xyz_content,
                "irc_reactant_xyz": irc_reactant_xyz_content,
                "irc_product_xyz": irc_product_xyz_content
            },
            "files": {
                "reactant_opt": str(work_dir / "reactant_opt.xyz"),
                "product_opt": str(work_dir / "product_opt.xyz"),
                "ts_opt": str(work_dir / "ts_opt.xyz"),
                "irc_reactant": str(result['irc_reactant_path']) if result.get('irc_reactant_path') else None,
                "irc_product": str(result['irc_product_path']) if result.get('irc_product_path') else None
            },
            "calculation_info": {
                "successful": result['success'],
                "error_message": result.get('error_message'),
                "timestamp": time.strftime('%Y-%m-%d %H:%M:%S')
            }
        }
        
        json_file = work_dir / "reaction_data.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
    
    def run_workflow(self):
        """运行完整的批量工作流"""
        logger.info("开始ORCA批量反应计算工作流")
        
        # 获取所有反应文件夹
        reaction_folders = self.get_reaction_folders()
        # 并行处理反应
        results = []
        skipped_count = 0
        
        if reaction_folders:
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                # 提交所有任务
                future_to_reaction = {
                    executor.submit(self.process_single_reaction, reaction_folder): reaction_folder
                    for reaction_folder in reaction_folders
                }
                
                # 收集结果
                for future in as_completed(future_to_reaction):
                    reaction_folder = future_to_reaction[future]
                    try:
                        result = future.result()
                        results.append(result)
                        
                        if result.get('skipped', False):
                            skipped_count += 1
                            if result['success']:
                                logger.info(f"⊙ {result['reaction_name']} 已完成（跳过）")
                            else:
                                logger.info(f"⊙ {result['reaction_name']} 跳过: {result.get('error_message', '未知原因')}")
                        elif result['success']:
                            logger.info(f"✓ {result['reaction_name']} 完成")
                        else:
                            logger.error(f"✗ {result['reaction_name']} 失败: {result.get('error_message', '未知错误')}")
                    except Exception as e:
                        logger.error(f"✗ {reaction_folder.name} 执行异常: {e}")

        if self.missing_inputs:
            results.extend(self.missing_inputs)

        if not results:
            logger.warning("未找到任何有效的反应文件夹")
            return
        
        # 生成总结报告
        self.generate_summary_report(results)
        self.save_results_csv(results)
        # 写入完成标志文件
        flag_name = f"done_{self.index}.flag" if self.index is not None else "done.flag"
        flag_path = self.output_base_dir / flag_name
        try:
            with open(flag_path, "w") as f:
                f.write(time.strftime('%Y-%m-%d %H:%M:%S'))
            logger.info(f"完成标志已写入 {flag_path}")
        except Exception as e:
            logger.warning(f"写入完成标志失败: {e}")
        
        logger.info(f"工作流完成，处理了 {len(results)} 个反应（跳过 {skipped_count} 个）")
    
    def generate_summary_report(self, results: List[Dict]):
        """生成总结报告"""
        report_file = self.output_base_dir / "batch_summary_report.txt"
        
        successful = [r for r in results if r['success'] and not r.get('skipped', False)]
        skipped = [r for r in results if r.get('skipped', False) and r['success']]
        failed = [r for r in results if not r['success']]
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("ORCA批量反应计算总结报告\n")
            f.write("=" * 50 + "\n")
            f.write(f"处理时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总反应数: {len(results)}\n")
            f.write(f"新计算成功数: {len(successful)}\n")
            f.write(f"跳过数（已完成）: {len(skipped)}\n")
            f.write(f"失败数: {len(failed)}\n")
            f.write("\n")
            
            if successful:
                f.write("新计算成功的反应:\n")
                f.write("-" * 50 + "\n")
                for result in successful:
                    f.write(f"反应: {result['reaction_name']}\n")
                    if result.get('original_folder'):
                        f.write(f"  原始文件夹: {result['original_folder']}\n")
                    if result['barrier_height']:
                        f.write(f"  能垒: {result['barrier_height']:.6f} Hartree ({result['barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    if result.get('reverse_barrier_height'):
                        f.write(f"  逆反应能垒: {result['reverse_barrier_height']:.6f} Hartree ({result['reverse_barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    if result['reaction_energy']:
                        f.write(f"  反应能: {result['reaction_energy']:.6f} Hartree ({result['reaction_energy'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    f.write(f"  虚频数: {sum(1 for freq in result['frequencies'] if freq < 0)}\n")
                    f.write(f"  有效过渡态: {'是' if result.get('is_valid_ts', False) else '否'}\n")
                    
                    # IRC分析信息
                    irc_analysis = result.get('irc_analysis', {})
                    if irc_analysis.get('performed', False):
                        f.write(f"  IRC计算: {'成功' if irc_analysis.get('success', False) else '失败'}\n")
                        if irc_analysis.get('success', False):
                            f.write(f"  反应物同构: {'是' if irc_analysis.get('reactant_isomorphic', False) else '否'}\n")
                            f.write(f"  产物同构: {'是' if irc_analysis.get('product_isomorphic', False) else '否'}\n")
                            f.write(f"  两端都同构: {'是' if irc_analysis.get('both_isomorphic', False) else '否'}\n")
                    else:
                        f.write(f"  IRC计算: 未执行（过渡态无效）\n")
                    f.write("\n")
            
            if skipped:
                f.write("跳过的反应（已完成）:\n")
                f.write("-" * 50 + "\n")
                for result in skipped:
                    f.write(f"反应: {result['reaction_name']}\n")
                    if result.get('original_folder'):
                        f.write(f"  原始文件夹: {result['original_folder']}\n")
                    if result.get('barrier_height'):
                        f.write(f"  能垒: {result['barrier_height']:.6f} Hartree ({result['barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    if result.get('reverse_barrier_height'):
                        f.write(f"  逆反应能垒: {result['reverse_barrier_height']:.6f} Hartree ({result['reverse_barrier_height'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    if result.get('reaction_energy'):
                        f.write(f"  反应能: {result['reaction_energy']:.6f} Hartree ({result['reaction_energy'] * HARTREE_TO_KCAL:.2f} kcal/mol)\n")
                    f.write("\n")
            
            if failed:
                f.write("失败的反应:\n")
                f.write("-" * 50 + "\n")
                for result in failed:
                    f.write(f"反应: {result['reaction_name']}\n")
                    if result.get('original_folder'):
                        f.write(f"  原始文件夹: {result['original_folder']}\n")
                    f.write(f"  错误: {result.get('error_message', '未知错误')}\n")
                    f.write("\n")
        
        # 生成JSON格式总结
        json_summary = {
            "summary": {
                "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
                "total_reactions": len(results),
                "successful_new": len(successful),
                "skipped": len(skipped),
                "failed": len(failed)
            },
            "reactions": []
        }
        
        for result in successful + skipped:
            reaction_info = {
                "reaction_name": result['reaction_name'],
                "original_folder": result.get('original_folder', ''),
                "barrier_height_hartree": result.get('barrier_height'),
                "barrier_height_kcal_mol": result['barrier_height'] * HARTREE_TO_KCAL if result.get('barrier_height') else None,
                "reverse_barrier_height_hartree": result.get('reverse_barrier_height'),
                "reverse_barrier_height_kcal_mol": result['reverse_barrier_height'] * HARTREE_TO_KCAL if result.get('reverse_barrier_height') else None,
                "reaction_energy_hartree": result.get('reaction_energy'),
                "reaction_energy_kcal_mol": result['reaction_energy'] * HARTREE_TO_KCAL if result.get('reaction_energy') else None,
                "imaginary_frequencies_count": sum(1 for freq in result.get('frequencies', []) if freq < 0),
                "is_valid_transition_state": result.get('is_valid_ts', False),
                "irc_analysis": result.get('irc_analysis', {
                    'performed': False,
                    'success': False,
                    'reactant_isomorphic': False,
                    'product_isomorphic': False,
                    'both_isomorphic': False
                }),
                "output_directory": str(result.get('work_dir', '')),
                "skipped": result.get('skipped', False)
            }
            json_summary["reactions"].append(reaction_info)
        
        json_file = self.output_base_dir / "batch_summary.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, ensure_ascii=False)
        
        logger.info(f"总结报告已保存: {report_file}")
        logger.info(f"JSON总结已保存: {json_file}")

    def save_results_csv(self, results: List[Dict[str, Any]]):
        """保存结果到csv"""
        if self.save_csv_path is None:
            self.save_csv_path = self.output_base_dir / "orca_validation_results.csv"
        self.save_csv_path.parent.mkdir(parents=True, exist_ok=True)

        fieldnames = [
            "row_index",
            "origin_idx",
            "origin_rsmi",
            "aug_idx",
            "working_dir",
            "save_dir",
            "reaction_name",
            "output_dir",
            "reactant_xyz",
            "product_xyz",
            "success",
            "skipped",
            "error_message",
            "reactant_energy_hartree",
            "product_energy_hartree",
            "ts_energy_hartree",
            "barrier_height_hartree",
            "reverse_barrier_height_hartree",
            "reaction_energy_hartree",
            "barrier_height_kcal_mol",
            "reverse_barrier_height_kcal_mol",
            "reaction_energy_kcal_mol",
            "imaginary_frequencies_count",
            "has_imaginary_freq",
            "is_valid_transition_state",
            "irc_performed",
            "irc_success",
            "irc_endpoint_match",
            "irc_rxn_status"
        ]

        sorted_results = sorted(
            results,
            key=lambda r: r.get("row_index") if r.get("row_index") is not None else -1
        )
        with open(self.save_csv_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for res in sorted_results:
                irc_analysis = res.get("irc_analysis") or {}
                barrier_height = res.get("barrier_height")
                reverse_barrier = res.get("reverse_barrier_height")
                reaction_energy = res.get("reaction_energy")
                row = {
                    "row_index": res.get("row_index"),
                    "origin_idx": res.get("origin_idx"),
                    "origin_rsmi": res.get("origin_rsmi"),
                    "aug_idx": res.get("aug_idx"),
                    "working_dir": res.get("working_dir"),
                    "save_dir": res.get("save_dir"),
                    "reaction_name": res.get("reaction_name"),
                    "output_dir": str(res.get("work_dir") or ""),
                    "reactant_xyz": res.get("reactant_src"),
                    "product_xyz": res.get("product_src"),
                    "success": res.get("success"),
                    "skipped": res.get("skipped"),
                    "error_message": res.get("error_message"),
                    "reactant_energy_hartree": res.get("reactant_energy"),
                    "product_energy_hartree": res.get("product_energy"),
                    "ts_energy_hartree": res.get("ts_energy"),
                    "barrier_height_hartree": barrier_height,
                    "reverse_barrier_height_hartree": reverse_barrier,
                    "reaction_energy_hartree": reaction_energy,
                    "barrier_height_kcal_mol": barrier_height * HARTREE_TO_KCAL if barrier_height is not None else None,
                    "reverse_barrier_height_kcal_mol": reverse_barrier * HARTREE_TO_KCAL if reverse_barrier is not None else None,
                    "reaction_energy_kcal_mol": reaction_energy * HARTREE_TO_KCAL if reaction_energy is not None else None,
                    "imaginary_frequencies_count": res.get("imaginary_freq_count", 0),
                    "has_imaginary_freq": res.get("has_imaginary_freq"),
                    "is_valid_transition_state": res.get("is_valid_ts"),
                    "irc_performed": irc_analysis.get("performed"),
                    "irc_success": irc_analysis.get("success"),
                    "irc_endpoint_match": irc_analysis.get("endpoint_match"),
                    "irc_rxn_status": irc_analysis.get("rxn_status")
                }
                writer.writerow(row)

        logger.info(f"结果csv已保存: {self.save_csv_path}")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='ORCA批量反应计算工作流（csv输入）')
    parser.add_argument('--csv', '-c', required=True,
                       help='反应对csv文件路径')
    parser.add_argument('--output', '-o', required=True,
                       help='输出基础目录')
    parser.add_argument('--orca', default="/inspire/hdd/project/chemicalreaction/misixuan-CZXS24220243/soft/orca_6_0_1_linux_x86-64_shared_openmpi416_avx2/orca",
                       help='ORCA可执行文件路径 (默认: orca)')
    parser.add_argument('--workers', '-w', type=int, default=1,
                       help='并行进程数 (默认: 1)')
    parser.add_argument('--force', '-f', action='store_true',
                       help='强制重新计算所有反应（忽略已有结果）')
    parser.add_argument('--save-csv', default=None, help='保存结果csv路径（默认: 输出目录/orca_validation_results.csv）')
    parser.add_argument('--index', type=int, default=None, help='仅处理csv中的指定行索引（从0开始）')
    parser.add_argument('--charge', type=int, default=0, help='默认电荷 (默认: 0)')
    parser.add_argument('--multiplicity', type=int, default=1, help='默认自旋多重度 (默认: 1)')
    
    args = parser.parse_args()
    
    try:
        # 创建工作流实例
        workflow = OrcaBatchWorkflow(
            csv_path=args.csv,
            output_base_dir=args.output,
            orca_exec=args.orca,
            max_workers=args.workers,
            force_recalculate=args.force,
            index=args.index,
            default_charge=args.charge,
            default_multiplicity=args.multiplicity,
            save_csv_path=args.save_csv
        )
        
        # 运行工作流
        workflow.run_workflow()
        
    except Exception as e:
        logger.error(f"工作流执行失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
