import pandas as pd
import os
import re  # 用于提取 `-` 后面的数字
import datetime

# 创建 Input 文件夹（如果不存在的话）
os.makedirs('./Input', exist_ok=True)

# 读取 CSV 数据
print("正在读取 CSV 数据...")
excel_data = pd.read_csv('./TNBC_immune_CODEX.csv')

# 提取需要的列：第一列（无列名）、patient, x, y, celltype, treatment
columns_needed = [excel_data.columns[0], 'patient', 'x', 'y', 'celltype', 'treatment']
data = excel_data[columns_needed]

print(f"数据集共包含 {len(data)} 行数据。")

# 提取 `-` 后的数字
def extract_number_after_dash(cell_string):
    match = re.search(r'-(\d+)', cell_string)  # 匹配 `-数字`
    return match.group(1) if match else "000000"  # 找不到返回默认值

# 解析 `h编号` 的排序方式
def parse_h_number(patient_id):
    h_part = patient_id.split('_')[0][7:]  # 获取 `h` 之后的部分（可能是 "03" 或 "03T2"）
    num_part = ''.join(filter(str.isdigit, h_part))  # 提取数字部分
    suffix = ''.join(filter(str.isalpha, h_part))  # 提取字母部分
    return (int(num_part), suffix)  # 返回 (数字部分, 字母后缀)

# **生成新的 patient_id**
new_patient_ids = []
print("\n正在生成 `patient_id`...")
for _, row in data.iterrows():
    h_number = row['patient'][1:]  # 直接去掉 "h"，保留 "03" 或 "03T2"
    treatment = int(row['treatment'])  # 读取 `treatment` 数值
    extracted_number = extract_number_after_dash(row[excel_data.columns[0]])  # 提取 `-` 后的数字

    new_patient_id = f"patient{h_number}_{treatment}_{extracted_number}"
    new_patient_ids.append(new_patient_id)

    print(f"  -> 提取 `patient_id`: {new_patient_id}")

# 添加 `new_patient_id` 列
data.loc[:, 'new_patient_id'] = new_patient_ids  # 避免 `SettingWithCopyWarning`

# **去重并按升序排序**
print("\n正在排序 `ImageNameList.txt`...")
unique_patient_ids = sorted(
    set(new_patient_ids),
    key=lambda x: (parse_h_number(x), int(x.split('_')[1]), int(x.split('_')[2]))
)

# **生成 ImageNameList.txt**
image_list_path = './Input/ImageNameList.txt'
with open(image_list_path, 'w') as image_name_list:
    for patient_id in unique_patient_ids:
        image_name_list.write(f"{patient_id}\n")
print(f"✅ `ImageNameList.txt` 已生成，共包含 {len(unique_patient_ids)} 个 `patient_id`。\n")

# **按 `patient_id` 存储数据**
print("正在创建 `Coordinates.txt`, `CellTypeLabel.txt`, `GraphLabel.txt` 文件...")
for idx, patient_id in enumerate(unique_patient_ids):
    # 过滤出当前 `patient_id` 的数据
    patient_data = data[data['new_patient_id'] == patient_id]

    # 获取坐标和细胞类型
    coordinates = patient_data[['x', 'y']].values
    celltypes = patient_data['celltype'].values
    treatment = patient_data['treatment'].values[0] - 1  # `treatment` 从 0 开始

    # 生成文件名
    coordinates_file = f'./Input/{patient_id}_Coordinates.txt'
    celltype_file = f'./Input/{patient_id}_CellTypeLabel.txt'
    graph_label_file = f'./Input/{patient_id}_GraphLabel.txt'

    # **保存坐标**
    with open(coordinates_file, 'w') as f:
        for coord in coordinates:
            f.write(f"{coord[0]}\t{coord[1]}\n")

    # **保存细胞类型**
    with open(celltype_file, 'w') as f:
        for celltype in celltypes:
            f.write(f"{celltype}\n")

    # **保存 treatment（图标签）**
    with open(graph_label_file, 'w') as f:
        f.write(f"{treatment}\n")

    print(f"  ({idx+1}/{len(unique_patient_ids)}) ✅ `{patient_id}` 数据已保存。")

# **结束**
print("\n🎉 所有病人的数据提取完成！")
print("📅 结束时间:", datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
