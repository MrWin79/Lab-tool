# -*- coding: utf-8 -*-
import os
import re
import io
import zipfile
import datetime
import gc 
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # 🛡️ 内存护盾：强制 Matplotlib 纯后台无头模式运行，防崩溃
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from openpyxl import load_workbook
import streamlit as st
import plotly.graph_objects as go
from streamlit_sortables import sort_items # 🖱️ 引入拖拽神级组件

# =====================================================
# Core Logic (Math)
# =====================================================
def parse_phases_from_wave(xls):
    for s in xls.sheet_names:
        if "assay" in s.lower():
            df = pd.read_excel(xls, sheet_name=s, header=None)
            target_aliases = {
                "Baseline": ["baseline"], "Etomoxir": ["eto"],
                "Oligo": ["oligo"], "FCCP": ["fccp"], "Rot": ["rot", "rra"]
            }
            phase_col_map = {}
            found_phases = []
            for i, row in df.iterrows():
                row_strs = [str(x).lower().strip() for x in row.values]
                if "baseline" in row_strs:
                    for col_idx, cell_val in enumerate(row_strs):
                        for std_name, aliases in target_aliases.items():
                            if std_name not in found_phases:
                                if any(a in cell_val for a in aliases):
                                    phase_col_map[std_name] = col_idx
                                    found_phases.append(std_name)
                                    break
                if phase_col_map and any("cycles" in str(x).lower() for x in row_strs):
                    parsed_any = False
                    result = {}
                    for p, col_idx in phase_col_map.items():
                        if col_idx < len(row.values):
                            cell_val = str(row.values[col_idx]).lower()
                            match = re.search(r'\d+', cell_val)
                            if match:
                                result[p] = int(match.group(0))
                                parsed_any = True
                    if parsed_any:
                        return [{"name": p, "cycles": result.get(p, 3)} for p in found_phases]
    return [{"name": "Baseline", "cycles": 3}, {"name": "Oligo", "cycles": 3},
            {"name": "FCCP", "cycles": 3}, {"name": "Rot", "cycles": 3}]

def derive_ranges_from_cycles(cycle_dict, phase_order):
    cur = 1
    ranges = {}
    for k in phase_order:
        if k in cycle_dict:
            ranges[k] = list(range(cur, cur + cycle_dict[k]))
            cur += cycle_dict[k]
    return ranges

def parse_cell_counts(xls):
    for s in xls.sheet_names:
        if "assay" in s.lower():
            df = pd.read_excel(xls, sheet_name=s, header=None)
            start_row, start_col = None, None
            for i, row in df.iterrows():
                row_strs = [str(x).lower().strip() for x in row.values]
                if "normalization values:" in row_strs:
                    start_row = i + 1  
                    for j, val in enumerate(df.iloc[i+1].values):
                        if str(val).strip().upper() == 'A':
                            start_col = j + 1  
                            break
                    break
            if start_row is None or start_col is None:
                start_row, start_col = 84, 2
            counts = {}
            rows_letters = "ABCDEFGH"
            for i, r_letter in enumerate(rows_letters):
                for j in range(12):
                    well = f"{r_letter}{j+1:02d}"
                    try:
                        val = df.iloc[start_row + i, start_col + j]
                        if pd.notna(val):
                            counts[well] = float(val)
                    except:
                        pass
            return counts
    return {}

def detect_wells_by_fg_theme(file_obj):
    wb = load_workbook(file_obj, data_only=True)
    sheet = next(s for s in wb.sheetnames if "assay" in s.lower())
    ws = wb[sheet]
    selected, unselected = set(), set()
    start_row, start_col = None, None
    for r in range(1, 30):
        for c in range(1, 10):
            cell_val = str(ws.cell(row=r, column=c).value).strip()
            if cell_val in ["A01", "A1"]:
                start_row, start_col = r, c
                break
        if start_row: break
    if not start_row: start_row, start_col = 12, 3 
    rows = "ABCDEFGH"
    for r_idx, r in enumerate(rows):
        for c_idx in range(12):
            cell = ws.cell(row=start_row + r_idx, column=start_col + c_idx)
            well = f"{r}{c_idx + 1:02d}"
            fg = cell.fill.fgColor if cell.fill else None
            if fg and fg.type == "theme" and fg.theme == 0:
                unselected.add(well)
            else:
                selected.add(well)
    return selected, unselected

def load_rate_sheet(xls, sheet_name, selected_wells):
    df = pd.read_excel(xls, sheet_name=sheet_name)
    df.columns = df.columns.str.lower()
    df = df[df["well"].isin(selected_wells)]
    df = df[~df["group"].isin(["Background", "Unassigned", "background", "unassigned"])]
    return df

def compute_metrics(df, ranges):
    df = df[~df["group"].isin(["Background", "Unassigned", "background", "unassigned"])].copy()
    last_base = max(ranges.get("Baseline", [1]))
    
    basal = df[df["measurement"] == last_base].groupby(["group", "well"], observed=True)["ocr"].mean()
    metrics = {}
    if "Oligo" in ranges:
        oligo = df[df["measurement"].isin(ranges["Oligo"])].groupby(["group", "well"], observed=True)["ocr"].min()
        metrics["ATP-linked Respiration"] = basal - oligo
    if "Rot" in ranges:
        rot = df[df["measurement"].isin(ranges["Rot"])].groupby(["group", "well"], observed=True)["ocr"].min()
        metrics["Basal Respiration"] = basal - rot
    if "Oligo" in ranges and "Rot" in ranges:
        metrics["Proton Leak"] = oligo - rot
    if "FCCP" in ranges and "Rot" in ranges:
        fccp = df[df["measurement"].isin(ranges["FCCP"])].groupby(["group", "well"], observed=True)["ocr"].max()
        metrics["Max Respiration"] = fccp - rot
    if "Etomoxir" in ranges:
        eto = df[df["measurement"].isin(ranges["Etomoxir"])].groupby(["group", "well"], observed=True)["ocr"].min()
        metrics["Fatty Acid Oxidation"] = basal - eto
        
    res = pd.DataFrame(metrics).reset_index().dropna(how='all', subset=list(metrics.keys()))
    res = res[~res["group"].isin(["Background", "Unassigned", "background", "unassigned"])]
    return res

def _row_letters_to_index(row_letters: str) -> int:
    idx = 0
    for ch in row_letters: idx = idx * 26 + (ord(ch) - ord('A') + 1)
    return idx

def _sanitize_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")


# =====================================================
# 🚀 极客优化：Zero-Disk I/O 纯内存写入引擎
# =====================================================
def save_plate_qc_to_zip(unselected, zipf):
    plate = np.zeros((8, 12))
    for r in range(8):
        for c in range(12):
            if f"{chr(ord('A') + r)}{c + 1:02d}" in unselected:
                plate[r, c] = 1
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.imshow(plate, cmap="gray_r")
    ax.set_xticks(range(12)); ax.set_xticklabels(range(1, 13))
    ax.set_yticks(range(8)); ax.set_yticklabels(list("ABCDEFGH"))
    ax.set_title("Plate QC (Grey = Excluded)")
    fig.tight_layout()
    
    img_bytes = io.BytesIO()
    fig.savefig(img_bytes, format='png', dpi=300)
    zipf.writestr("plate_qc.png", img_bytes.getvalue())
    plt.close(fig)

def export_prism_csvs_to_zip(summary_df, zipf, base_dir, suffix):
    df = summary_df.copy()
    df.columns = [c.strip() for c in df.columns]
    metric_cols = [c for c in df.columns if c not in ["group", "well"]]
    if not metric_cols: return
    row_letters = df["well"].astype(str).str.extract(r"^([A-Za-z]+)")[0].str.upper()
    col_nums = df["well"].astype(str).str.extract(r"(\d+)$")[0]
    df["_row_idx"] = row_letters.apply(_row_letters_to_index)
    df["_col_num"] = col_nums.astype(int)
    
    for metric in metric_cols:
        tmp = df[["group", "well", "_row_idx", "_col_num", metric]].copy()
        tmp = tmp.sort_values(["group", "_row_idx", "_col_num", "well"], kind="mergesort")
        tmp["_rep"] = tmp.groupby("group", observed=True).cumcount() + 1
        wide = tmp.pivot(index="_rep", columns="group", values=metric).sort_index()
        wide = wide.rename_axis(None, axis=0).rename_axis(None, axis=1)
        csv_str = wide.to_csv(index=False)
        # 写入带有BOM的UTF8，防止Excel乱码
        zipf.writestr(f"{base_dir}/prism_csv/Prism_{_sanitize_filename(metric)}_{suffix}.csv", csv_str.encode("utf-8-sig"))

def save_results_to_zip(df, zipf, base_dir, suffix, control):
    csv_str = df.to_csv(index=False)
    zipf.writestr(f"{base_dir}/summary_{suffix}.csv", csv_str.encode("utf-8-sig"))
    export_prism_csvs_to_zip(df, zipf, base_dir, suffix)

    for metric in df.columns[2:]:
        fig, ax = plt.subplots(figsize=(8, 5))
        g = df.groupby("group", observed=True)[metric]
        m, s = g.mean(), g.std()
        ax.bar(m.index, m.values, yerr=s.values, capsize=4, zorder=1, alpha=0.8, color="#4A90E2")
        ax.set_title(metric, fontweight="bold")

        for i, grp in enumerate(m.index):
            y_vals = df[df["group"] == grp][metric].dropna()
            if len(y_vals) > 0:
                x_vals = np.random.normal(i, 0.06, size=len(y_vals))
                ax.scatter(x_vals, y_vals, color='black', alpha=0.6, s=25, zorder=3)

        ylabel = "Cell Count" if metric == "Cell Number" else "OCR (pmol/min)" if suffix == "raw" else "OCR (pmol/min/Norm. Unit)"
        ax.set_ylabel(ylabel)

        if control in m.index:
            ctrl_vals = df[df["group"] == control][metric].dropna()
            valid_maxes = [df[metric].max()]
            if not m.isna().all() and not s.isna().all(): valid_maxes.append((m + s).max())
            global_ymax = max([v for v in valid_maxes if not np.isnan(v)] or [0])
            for i, grp in enumerate(m.index):
                if grp == control: continue
                vals = df[df["group"] == grp][metric].dropna()
                if len(ctrl_vals) >= 2 and len(vals) >= 2:
                    _, p = ttest_ind(ctrl_vals, vals, equal_var=False)
                    stars = "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    p_str = "p<0.0001" if p < 0.0001 else f"p={round(p, 4)}"
                    grp_max_val = max(m[grp] + s[grp], vals.max())
                    ax.text(i, grp_max_val + global_ymax * 0.03, f"{stars}\n{p_str}", ha="center", va="bottom", fontsize=8, linespacing=1.2)

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        y_min = df[metric].min()
        if not np.isnan(y_min):
            ax.set_ylim(min(0, y_min * 1.1) if y_min < 0 else 0, ax.get_ylim()[1] * 1.25)

        fig.tight_layout()
        img_bytes = io.BytesIO()
        fig.savefig(img_bytes, format='png', dpi=300)
        zipf.writestr(f"{base_dir}/{metric}_{suffix}.png", img_bytes.getvalue())
        plt.close(fig)
        
    plt.close('all')

def save_kinetic_graphs_to_zip(df, zipf, base_dir, suffix, ranges, phase_order, metric_col="ocr"):
    if metric_col not in df.columns: return
    g = df.groupby(["group", "measurement"], observed=True)[metric_col]
    m, s = g.mean().unstack(level=0), g.sem().unstack(level=0)
    if m.empty: return

    fig, ax = plt.subplots(figsize=(10, 5)) 
    for col in m.columns:
        ax.errorbar(m.index, m[col], yerr=s[col], label=col, marker='o', capsize=3, markersize=5, linewidth=2)

    ymin, ymax = ax.get_ylim()
    for i in range(1, len(phase_order)):
        phase = phase_order[i]
        if phase in ranges and len(ranges[phase]) > 0:
            line_x = ranges[phase][0] - 0.5
            ax.axvline(x=line_x, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
            ax.text(line_x, ymax * 0.95, phase, rotation=0, ha='center', va='top', fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))

    metric_upper = metric_col.upper()
    ax.set_title(f"Kinetic {metric_upper} Profile", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("Measurement", fontsize=12
