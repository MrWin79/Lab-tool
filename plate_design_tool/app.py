import streamlit as st
import pandas as pd
import random
import time
import math
import io
from collections import defaultdict

# ==========================================
# 🎨 页面配置与 UI 样式
# ==========================================
st.set_page_config(page_title="TDIMP Plate Designer | 智能微孔板阵列系统", layout="wide", page_icon="🧪")

COLORS = [
    "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", 
    "#E6B3FF", "#FFB3E6", "#E0E0E0", "#FFC8A2", "#D4F0F0",
    "#FF9CEE", "#C5A3FF", "#BFFCC6", "#FFC9DE", "#D5AAFF"
]

CTRL_COLORS = [
    "#1e3a8a", "#047857", "#b45309", "#be123c", "#4338ca", 
    "#0f766e", "#6d28d9", "#a21caf", "#1d4ed8", "#15803d"
]

def apply_color(val, gene_keys, ctrl_keys):
    if pd.isna(val) or val == "": return "" 
    if "[排除]" in str(val): return "color: #64748b;" 
    if "⚠️" in str(val): return "background-color: #7f1d1d; color: #ffffff;"
        
    clean_val = str(val).strip()
    if clean_val in ctrl_keys:
        idx = ctrl_keys.index(clean_val) % len(CTRL_COLORS)
        return f"background-color: {CTRL_COLORS[idx]}; color: #ffffff;"
    if clean_val in gene_keys:
        idx = gene_keys.index(clean_val) % len(COLORS)
        return f"background-color: {COLORS[idx]}; color: #000000;"
    return ""

def to_excel(df):
    output = io.BytesIO()
    # 完美兼容 Pandas 2.1+
    clean_df = df.map(lambda x: str(x).replace(" ⚠️", "") if isinstance(x, str) else x)
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        clean_df.to_excel(writer, index=True, sheet_name='PlateMap')
    return output.getvalue()

# ==========================================
# 🧠 第一代核心：10000次蒙特卡洛纯净寻优引擎 (Monte Carlo)
# ==========================================
def get_quadrant(r_idx, c_idx, total_rows, total_cols):
    is_bottom = r_idx >= total_rows // 2
    is_right = c_idx >= total_cols // 2
    if not is_bottom and not is_right: return 1
    if not is_bottom and is_right: return 2
    if is_bottom and not is_right: return 3
    return 4

def get_sq_dist(r1, c1, r2, c2):
    return (ord(r1) - ord(r2))**2 + (int(c1) - int(c2))**2

def is_quad_valid(quad, q_counts, req_reps):
    if req_reps == 0: return False
    max_per_quad = math.ceil(req_reps / 4.0)
    num_max_quads = req_reps % 4 if req_reps % 4 != 0 else 4
    proposed_count = q_counts.get(quad, 0) + 1
    if proposed_count > max_per_quad: return False
    if proposed_count == max_per_quad:
        if sum(1 for v in q_counts.values() if v == max_per_quad) >= num_max_quads:
            return False
    return True

def run_simulation(genes_dict, ctrls_dict, all_rows, all_cols, available_wells, seed):
    random.seed(seed)
    
    layout_result = {}
    gene_placed_coords = defaultdict(list)
    ctrl_placed_coords = defaultdict(list)
    
    gene_keys = list(genes_dict.keys())
    random.shuffle(gene_keys)
    
    MAX_ATTEMPTS = 10000
    t_start = time.perf_counter()
    success = False
    strategy = ""
    attempt = 0
    
    for attempt in range(1, MAX_ATTEMPTS + 1):
        # 阶梯式降级策略阈值
        is_distance_greedy = (attempt <= 3000)
        require_gene_no_corner = (attempt <= 6000)
        require_ctrl_no_corner = (attempt <= 6000)
        require_ctrl_no_edge = (attempt <= 8000)
        require_ctrl_quad = (attempt <= 9000)
        
        available = set(available_wells)
        layout_result.clear()
        
        gene_placed_rows = defaultdict(set)
        gene_placed_cols = defaultdict(set)
        gene_placed_coords.clear()
        ctrl_placed_coords.clear()
        gene_placed_quads = defaultdict(lambda: defaultdict(int))
        ctrl_placed_quads = defaultdict(lambda: defaultdict(int))
        
        try:
            # --- 1. 靶点排布 (同行同列隔离绝对不可被降级) ---
            for g_name in gene_keys:
                req_reps = int(genes_dict[g_name])
                for _ in range(req_reps):
                    cands = []
                    for (r, c) in available:
                        if r not in gene_placed_rows[g_name] and c not in gene_placed_cols[g_name]:
                            # 靶点防共角
                            if require_gene_no_corner:
                                if any(get_sq_dist(r, c, pr, pc) <= 2 for pr, pc in gene_placed_coords[g_name]):
                                    continue
                            
                            quad = get_quadrant(all_rows.index(r), all_cols.index(c), len(all_rows), len(all_cols))
                            if is_quad_valid(quad, gene_placed_quads[g_name], req_reps):
                                d_sq = float('inf')
                                if gene_placed_coords[g_name]:
                                    d_sq = min(get_sq_dist(r, c, pr, pc) for pr, pc in gene_placed_coords[g_name])
                                cands.append((r, c, d_sq, quad))
                    
                    if not cands: raise ValueError("Gene Blocked")
                    
                    if is_distance_greedy:
                        max_d = max(cand[2] for cand in cands)
                        best = [cand for cand in cands if cand[2] == max_d]
                        chosen = random.choice(best)
                    else:
                        chosen = random.choice(cands)
                        
                    r, c, _, quad = chosen
                    available.remove((r, c))
                    layout_result[(r, c)] = g_name
                    gene_placed_rows[g_name].add(r)
                    gene_placed_cols[g_name].add(c)
                    gene_placed_coords[g_name].append((r, c))
                    gene_placed_quads[g_name][quad] += 1
            
            # --- 2. 对照组排布 ---
            for ctrl_name, req_reps in ctrls_dict.items():
                req_reps = int(req_reps)
                for _ in range(req_reps):
                    cands = []
                    for (r, c) in available:
                        quad = get_quadrant(all_rows.index(r), all_cols.index(c), len(all_rows), len(all_cols))
                        
                        if require_ctrl_quad and not is_quad_valid(quad, ctrl_placed_quads[ctrl_name], req_reps):
                            continue
                        
                        d_sq_same = float('inf')
                        if ctrl_placed_coords[ctrl_name]:
                            d_sq_same = min(get_sq_dist(r, c, pr, pc) for pr, pc in ctrl_placed_coords[ctrl_name])
                        
                        # 同类对照组一票否决制防共角/共边
                        if require_ctrl_no_corner and d_sq_same <= 2:
                            continue
                        if require_ctrl_no_edge and d_sq_same <= 1:
                            continue
                            
                        cands.append((r, c, d_sq_same, quad)) 
                    
                    if not cands: raise ValueError("Ctrl Blocked")
                    
                    if is_distance_greedy:
                        max_d = max(cand[2] for cand in cands)
                        best = [cand for cand in cands if cand[2] == max_d]
                        chosen = random.choice(best)
                    else:
                        chosen = random.choice(cands)
                    
                    r, c, _, quad = chosen
                    available.remove((r, c))
                    layout_result[(r, c)] = ctrl_name
                    ctrl_placed_coords[ctrl_name].append((r, c))
                    ctrl_placed_quads[ctrl_name][quad] += 1
            
            success = True
            
            # 记录当前满足的最高策略维度
            strat_list = []
            if is_distance_greedy: strat_list.append("最大空间拉扯")
            if require_gene_no_corner: strat_list.append("靶点不共角")
            else: strat_list.append("靶点降级(共角)")
            if require_ctrl_no_corner: strat_list.append("同类对照不共角")
            elif require_ctrl_no_edge: strat_list.append("同类对照不共边")
            else: strat_list.append("对照降级(紧凑)")
            strategy = " | ".join(strat_list)
            break
            
        except ValueError:
            continue

    elapsed = time.perf_counter() - t_start
    return success, layout_result, gene_placed_coords, ctrl_placed_coords, attempt, elapsed, strategy

# ==========================================
# 🖥️ Streamlit 界面逻辑
# ==========================================
st.title("🧪 TDIMP 核心阵列系统 | 蒙特卡洛纯净引擎")

if "plate_type" not in st.session_state: st.session_state.plate_type = "96孔板"
if "seed_val" not in st.session_state: st.session_state.seed_val = 42

with st.sidebar:
    st.header("🧫 板型设置")
    new_plate_type = st.radio("微孔板规格", ["96孔板", "384孔板"], index=0 if st.session_state.plate_type == "96孔板" else 1)
    if new_plate_type != st.session_state.plate_type:
        st.session_state.plate_type = new_plate_type
        if "plate_mask" in st.session_state: del st.session_state.plate_mask
        st.rerun()

is_96 = (st.session_state.plate_type == "96孔板")
all_rows = [chr(65+i) for i in range(8 if is_96 else 16)]
all_cols = [f"{i:02d}" for i in range(1, 13 if is_96 else 25)]

# 智能初始化默认空间掩码
if "plate_mask" not in st.session_state:
    mask = pd.DataFrame(True, index=all_rows, columns=all_cols)
    if is_96:
        # 96孔板默认排外1圈
        mask.iloc[0, :] = mask.iloc[-1, :] = mask.iloc[:, 0] = mask.iloc[:, -1] = False
    else:
        # 384孔板默认排外2圈
        mask.iloc[0:2, :] = mask.iloc[-2:, :] = mask.iloc[:, 0:2] = mask.iloc[:, -2:] = False
    st.session_state.plate_mask = mask

tab1, tab2 = st.tabs(["🧬 阵列设计与排布运算", "🖱️ 可用空间精细控制 (Mask编辑器)"])

with tab2:
    st.subheader("可用物理空间编辑器")
    st.caption("提示：可以通过上方【批量快捷操作】快速圈定整行/整列，再在下方的网格图表中进行单击微调。")
    
    # ---------------- ⚡ 批量快捷操作面板 ----------------
    with st.container():
        st.markdown("#### ⚡ 批量快捷操作")
        col_r, col_c, col_btn = st.columns([2, 2, 1.5])
        
        with col_r:
            sel_rows = st.multiselect("框选目标【行】 (如 A, B, O, P)", all_rows)
        with col_c:
            sel_cols = st.multiselect("框选目标【列】 (如 01, 02, 23)", all_cols)
        with col_btn:
            st.write("") # 占位对齐
            c_b1, c_b2 = st.columns(2)
            with c_b1:
                if st.button("✔️ 设为可用", use_container_width=True):
                    if sel_rows: st.session_state.plate_mask.loc[sel_rows, :] = True
                    if sel_cols: st.session_state.plate_mask.loc[:, sel_cols] = True
                    st.rerun()
            with c_b2:
                if st.button("❌ 设为排除", use_container_width=True):
                    if sel_rows: st.session_state.plate_mask.loc[sel_rows, :] = False
                    if sel_cols: st.session_state.plate_mask.loc[:, sel_cols] = False
                    st.rerun()
                    
    st.divider()
    edited_mask = st.data_editor(st.session_state.plate_mask, use_container_width=True, height=350 if is_96 else 600)
    st.session_state.plate_mask = edited_mask

available_wells = [(r, c) for r in all_rows for c in all_cols if edited_mask.at[r, c]]
total_active_wells = len(available_wells)

with tab1:
    if "df_genes" not in st.session_state:
        st.session_state.df_genes = pd.DataFrame({"靶点名称 (Gene)": [f"Gene_{i}" for i in range(1, 11)], "复孔数量": [5] * 10})
    if "df_ctrls" not in st.session_state:
        st.session_state.df_ctrls = pd.DataFrame({"对照名称 (Control)": ["NTC"], "复孔数量": [6]})

    col_L, col_R = st.columns(2)
    with col_L:
        st.subheader("1. 实验靶点 (Genes)")
        edited_genes = st.data_editor(st.session_state.df_genes, num_rows="dynamic", use_container_width=True)
    with col_R:
        st.subheader("2. 对照组 (Controls)")
        edited_ctrls = st.data_editor(st.session_state.df_ctrls, num_rows="dynamic", use_container_width=True)

    c_info, c_seed, c_btn1, c_btn2 = st.columns([1.5, 1, 1.5, 1.5])
    with c_info: st.info(f"💡 当前可用空间: **{total_active_wells}** 孔")
    with c_seed: manual_seed = st.number_input("🎲 随机种子", value=st.session_state.seed_val, step=1, label_visibility="collapsed")
    st.session_state.seed_val = manual_seed
    
    run_algorithm = False
    with c_btn1:
        if st.button("🚀 启动蒙特卡洛寻优", type="primary", use_container_width=True): run_algorithm = True
    with c_btn2:
        if st.button("🔄 一键重新演化", use_container_width=True):
            st.session_state.seed_val = random.randint(1, 99999)
            run_algorithm = True

    if run_algorithm:
        genes_dict = {row["靶点名称 (Gene)"]: int(row["复孔数量"]) for _, row in edited_genes.iterrows() if pd.notna(row["靶点名称 (Gene)"])}
        ctrls_dict = {row["对照名称 (Control)"]: int(row["复孔数量"]) for _, row in edited_ctrls.iterrows() if pd.notna(row["对照名称 (Control)"])}
        total_wells_req = sum(genes_dict.values()) + sum(ctrls_dict.values())
        
        if total_wells_req > total_active_wells:
            st.error(f"❌ 空间严重不足！分配 {total_wells_req} 孔，但可用仅 {total_active_wells} 孔。")
        else:
            with st.spinner("🧠 万次级蒙特卡洛碰撞计算中..."):
                success, layout, gene_coords, ctrl_coords, attempts, elapsed, strategy = run_simulation(
                    genes_dict, ctrls_dict, all_rows, all_cols, available_wells, seed=st.session_state.seed_val
                )
                
            if not success:
                st.error("❌ 10000 次碰撞后未能求解。排布空间过度锁死，建议增加物理可用区域，或减少复孔数量。")
            else:
                conflicts = set()
                err_row_col, err_gene_corner, err_gene_quad, err_ctrl_quad, err_ctrl_corner = [], [], [], [], []

                for g, coords in gene_coords.items():
                    for r, c in coords:
                        same_r = [wc for wr, wc in coords if wr == r]
                        same_c = [wr for wr, wc in coords if wc == c]
                        if len(same_r) > 1: 
                            err_row_col.append(f"**{g}** 行冲突 ({r}行: {', '.join(same_r)})")
                            conflicts.update([(r, xc) for xc in same_r])
                        if len(same_c) > 1: 
                            err_row_col.append(f"**{g}** 列冲突 ({c}列: {', '.join(same_c)})")
                            conflicts.update([(xr, c) for xr in same_c])
                err_row_col = list(set(err_row_col))

                for g, coords in gene_coords.items():
                    for i, (r1, c1) in enumerate(coords):
                        for r2, c2 in coords[i+1:]:
                            r1_i, c1_i = all_rows.index(r1), all_cols.index(c1)
                            r2_i, c2_i = all_rows.index(r2), all_cols.index(c2)
                            if (r1_i - r2_i)**2 + (c1_i - c2_i)**2 <= 2:
                                err_gene_corner.append(f"**{g}** 局部接触 ({r1}{c1} ↔ {r2}{c2})")
                                conflicts.update([(r1, c1), (r2, c2)])

                for g, coords in gene_coords.items():
                    quads = [get_quadrant(all_rows.index(r), all_cols.index(c), len(all_rows), len(all_cols)) for r, c in coords]
                    q_counts = [quads.count(q) for q in [1,2,3,4]]
                    if max(q_counts) - min(q_counts) > 1:
                        err_gene_quad.append(f"**{g}** 分布失衡 (Q1:{q_counts[0]} Q2:{q_counts[1]} Q3:{q_counts[2]} Q4:{q_counts[3]})")

                for c_name, coords in ctrl_coords.items():
                    quads = [get_quadrant(all_rows.index(r), all_cols.index(c), len(all_rows), len(all_cols)) for r, c in coords]
                    q_counts = [quads.count(q) for q in [1,2,3,4]]
                    if max(q_counts) - min(q_counts) > 1:
                        err_ctrl_quad.append(f"**{c_name}** 分布失衡")

                for c_name, coords in ctrl_coords.items():
                    for i, (r1, c1) in enumerate(coords):
                        for r2, c2 in coords[i+1:]:
                            r1_i, c1_i = all_rows.index(r1), all_cols.index(c1)
                            r2_i, c2_i = all_rows.index(r2), all_cols.index(c2)
                            if (r1_i - r2_i)**2 + (c1_i - c2_i)**2 <= 2:
                                err_ctrl_corner.append(f"**{c_name}** 内部拥挤 ({r1}{c1} ↔ {r2}{c2})")
                                conflicts.update([(r1, c1), (r2, c2)])

                # ---------------- 🏆 核心加权打分系统 ----------------
                base_score = 100
                deduct_row_col = len(err_row_col) * 30
                deduct_gene_quad = len(err_gene_quad) * 20
                deduct_gene_corner = len(err_gene_corner) * 15
                deduct_ctrl_quad = len(err_ctrl_quad) * 10
                deduct_ctrl_corner = len(err_ctrl_corner) * 15
                
                final_score = base_score - (deduct_row_col + deduct_gene_quad + deduct_gene_corner + deduct_ctrl_quad + deduct_ctrl_corner)
                final_score = max(0, final_score)
                
                if final_score == 100: grade, grade_color = "S 级 (完美无瑕)", "#10b981"
                elif final_score >= 85: grade, grade_color = "A 级 (轻微妥协)", "#3b82f6"
                elif final_score >= 70: grade, grade_color = "B 级 (勉强可用)", "#f59e0b"
                else: grade, grade_color = "F 级 (建议重排)", "#ef4444"

                # ---------------- 构建静态 Pandas 矩阵 ----------------
                plate_matrix = []
                for r in all_rows:
                    row_data = []
                    for c in all_cols:
                        if (r, c) not in available_wells:
                            row_data.append("[排除]")
                        else:
                            val = layout.get((r, c), "")
                            if (r, c) in conflicts and val != "":
                                val += " ⚠️"
                            row_data.append(val)
                    plate_matrix.append(row_data)
                    
                df_plate = pd.DataFrame(plate_matrix, index=all_rows, columns=all_cols)
                gene_keys = list(gene_coords.keys())
                ctrl_keys = list(ctrl_coords.keys())
                
                cell_width = "40px" if not is_96 else "65px"
                cell_font = "10px" if not is_96 else "14px"

                st.success(f"✅ 排布计算完成！耗时 {elapsed:.3f} 秒，碰撞迭代 {attempts} 次。")
                
                st.markdown(f"""
                <div style="background-color: #f8fafc; border-radius: 10px; padding: 20px; border-left: 8px solid {grade_color}; margin-bottom: 20px;">
                    <h3 style="margin: 0; color: #1e293b;">🧬 阵列综合评分: <span style="color: {grade_color}; font-size: 1.5em;">{final_score}</span> / 100</h3>
                    <p style="margin: 5px 0 0 0; color: #64748b; font-weight: bold;">评级: {grade} | 当前收敛策略: {strategy}</p>
                </div>
                """, unsafe_allow_html=True)
                
                styled_plate = df_plate.style.map(
                    lambda x: apply_color(x, gene_keys, ctrl_keys)
                ).set_properties(**{
                    'text-align': 'center', 'border': '1px solid #e2e8f0',
                    'min-width': cell_width, 'font-size': cell_font
                })
                
                if conflicts:
                    st.warning("🚨 注意：当前满载状态触发了约束降级，图中标红并带有 ⚠️ 的孔位为冲突发生地。")
                
                st.markdown(f"### 🔬 {st.session_state.plate_type} 阵列全景视图")
                st.dataframe(styled_plate, use_container_width=True)

                col_dl, _ = st.columns([1, 2])
                with col_dl:
                    st.download_button(
                        label="📥 下载当前板图 (Excel)",
                        data=to_excel(df_plate),
                        file_name=f"TDIMP_PlateMap_Score{final_score}_{int(time.time())}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True
                    )
                
                st.markdown("---")
                st.markdown("### 📊 排布质量五维扫描诊断")
                c1, c2, c3, c4, c5 = st.columns(5)
                
                def render_metric(col, title, pass_text, err_list, delta_color):
                    if not err_list:
                        col.metric(title, "✅ 通过", pass_text)
                    else:
                        col.metric(title, "⚠️ 扣分项", f"{len(err_list)} 处异常", delta_color=delta_color)
                        with col.expander("🔍 详情"):
                            st.markdown("\n".join([f"- {e}" for e in err_list]))

                render_metric(c1, "1. 行列绝对隔离 (-30/个)", "0 冲突", err_row_col, "inverse")
                render_metric(c2, "2. 靶点防共角 (-15/个)", "0 接触", err_gene_corner, "inverse")
                render_metric(c3, "3. 靶点象限平滑 (-20/个)", "平滑分布", err_gene_quad, "inverse")
                render_metric(c4, "4. 各对照组平滑 (-10/个)", "平滑分布", err_ctrl_quad, "inverse")
                render_metric(c5, "5. 同类对照防聚集 (-15/个)", "0 扎堆", err_ctrl_corner, "off")

                st.markdown("<div style='margin-top: 30px; text-align: right; color: #94a3b8; font-size: 0.8em; font-family: monospace;'>Powered by TDIMP Monte Carlo Engine | Designed by Fxq</div>", unsafe_allow_html=True)
