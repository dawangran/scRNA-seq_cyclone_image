import os
import json
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from rapidfuzz.distance import Levenshtein
import argparse
import glob

# 按基因个数排序条形码
def sort_sequences_by_gene_count(data):
    """根据基因对应的条形码数量对基因进行降序排序"""
    return tuple(data.items())

# 提取基因和条形码个数，并计数
def extract_genes(data):
    """统计每个基因对应的条形码数量"""
    return defaultdict(int, {gene: len(sequences) for gene, sequences in data.items()})

# 按条形码个数排序基因
def sort_genes_by_sequence_count(gene_counts):
    """根据条形码数量对基因排序"""
    return sorted(gene_counts.items(), key=lambda item: item[1])

# 计算编辑距离（优化后）
def calculate_levenshtein_distance(barcode, s_barcode_list, max_distance):
    """计算与目标条形码在设定阈值范围内的编辑距离"""
    return [(s_barcode, Levenshtein.distance(s_barcode, barcode))
            for s_barcode in s_barcode_list
            if Levenshtein.distance(s_barcode, barcode) <= max_distance]



def calculate_result(gene_dict):
    # 初始化分类变量
    invalid_barcode=[]
    consistent_barcode = None
    wrong_barcode = {}  # 这个变量在原始代码中没有使用，如果需要，请添加相应逻辑
    multiple_barcode = {}
    best_barcode_result = []

    # 统计长度为1的值的出现频率
    length_one_count = {}

    # 遍历字典，分类键值对并统计长度为1的值的出现频率
    for key, value in gene_dict.items():
        if len(value) > 1:
            multiple_barcode[key] = value
        elif len(value) == 0:
            invalid_barcode.append(key)
        elif len(value) == 1:
            value_tuple = tuple(value)  # 将列表转换为元组，以便可以用作字典的键
            if value_tuple in length_one_count:
                length_one_count[value_tuple].append(key)
            else:
                length_one_count[value_tuple] = [key]

    # 找出出现频率最高的值
    if length_one_count:
        most_frequent_value = max(length_one_count, key=lambda k: len(length_one_count[k]))
        most_frequent_keys = length_one_count[most_frequent_value]
        consistent_barcode_len = len(most_frequent_keys)

        # 设置出现频率最高的值
        consistent_barcode = most_frequent_value[0]

        # 构建其余长度为1的值的字典
        for value, keys in length_one_count.items():
            if value != most_frequent_value:
                for key in keys:
                    wrong_barcode[key] = list(value)

    # 计算最佳条形码结果
    if consistent_barcode is not None:
        best_barcode_result = [consistent_barcode,consistent_barcode_len,len(wrong_barcode),len(multiple_barcode),len(invalid_barcode),len(gene_dict)]
        

    return {
        "invalid_barcode": invalid_barcode,
        "consistent_barcode": consistent_barcode,
        "wrong_barcode": wrong_barcode,  # 如果需要，请添加相应逻辑
        "multiple_barcode": multiple_barcode,
        "best_barcode_result": best_barcode_result
    }        


def best_barcode(barcode_gene, gene_barcode_data, min_dis):
    l_barcode, l_gene_list = barcode_gene
    gene_dict = {}
    for l_gene in l_gene_list:
        s_barcode_list = gene_barcode_data.get(l_gene, [])
        if s_barcode_list:
            if l_barcode in s_barcode_list:
                gene_dict[l_gene] = [l_barcode]
            else:
                dis_list = calculate_levenshtein_distance(l_barcode, s_barcode_list, min_dis)
                if dis_list:
                    min_dis_pair = min(dis_list, key=lambda x: x[1])
                    gene_dict[l_gene]=[t[0] for t in dis_list if t[1] == min_dis_pair[1]]
                else:
                    gene_dict[l_gene]=[]   
        else:
            continue
    return {l_barcode:calculate_result(gene_dict)}


# 批量处理条形码
def process_batch(batch, gene_barcode_data, min_dis):
    """并行处理条形码批次"""
    batch_results = {}
    for barcode_gene in batch:
        result = best_barcode(barcode_gene, gene_barcode_data, min_dis)
        batch_results.update(result)
    return batch_results


def merge_json_files(output_dir, merged_file_path,sep="best_dict"):
    merged_data = {}
    # 构建匹配所有批次文件的模式
    batch_files_pattern = os.path.join(output_dir, "barcode_{}_batch_*.json".format(sep))
    # 使用glob模块找到所有匹配的文件
    batch_files = glob.glob(batch_files_pattern)
    
    for file_path in batch_files:
        # 读取 JSON 文件内容
        with open(file_path, 'r') as f:
            data = json.load(f)
        # 合并数据
        merged_data.update(data)
    
    # 将合并后的数据保存到新的 JSON 文件中
    with open(merged_file_path, 'w') as f:
        json.dump(merged_data, f, indent=4)
    
    # 删除所有批次文件
    for file_path in batch_files:
        os.remove(file_path)

def process_barcodes(result, s_gb_dict, l_gb_dict):
    multi_dict = {}
    wrong_dict = {}
    best_barcode_dict = {}

    for k, v in result.items():
        best_barcode_list = v.get('best_barcode_result', [])
        lb = k
        if best_barcode_list:
            if v.get('multiple_barcode'):
                multiple_barcode = v['multiple_barcode']
                process_multiple_barcode(multiple_barcode, lb, best_barcode_list, multi_dict, s_gb_dict, l_gb_dict)

            if v.get('wrong_barcode'):
                wrong_barcode = v['wrong_barcode']
                process_wrong_barcode(wrong_barcode, lb, best_barcode_list, wrong_dict, s_gb_dict, l_gb_dict)

            score = best_barcode_list[1] / best_barcode_list[-1] if best_barcode_list[-1] != 0 else 0
            if score >= 0.6:
                best_barcode_dict[lb] = best_barcode_list[0]
        else:
            if v.get('multiple_barcode'):
                multiple_barcode = v['multiple_barcode']
                process_no_best_barcode(multiple_barcode, lb, best_barcode_dict, s_gb_dict, l_gb_dict)

    return multi_dict, wrong_dict, best_barcode_dict

def process_multiple_barcode(multiple_barcode, lb, best_barcode_list, multi_dict, s_gb_dict, l_gb_dict):
    for g, sb_list in multiple_barcode.items():
        lbu = f"{g}_{lb}"
        tmp_result = {}
        for k in sb_list:
            sbu = f"{g}_{k}"
            if sbu in s_gb_dict:
                u_dis = 10
                for lu in l_gb_dict.get(lbu, []):
                    for su in s_gb_dict.get(sbu, []):
                        if Levenshtein.distance(lu, su) <= u_dis:
                            u_dis = Levenshtein.distance(lu, su)
                if u_dis <= 4:
                    tmp_result[k] = u_dis
        if tmp_result:
            min_pair = min(tmp_result, key=tmp_result.get)
            update_best_barcode_list(best_barcode_list, min_pair, multi_dict, lb, g)

def process_wrong_barcode(wrong_barcode, lb, best_barcode_list, wrong_dict, s_gb_dict, l_gb_dict):
    for g, sb_list in wrong_barcode.items():
        sb_list.append(lb)
        lbu = f"{g}_{lb}"
        tmp_result = {}
        for k in sb_list:
            sbu = f"{g}_{k}"
            if sbu in s_gb_dict:
                u_dis = 10
                for lu in l_gb_dict.get(lbu, []):
                    for su in s_gb_dict.get(sbu, []):
                        if Levenshtein.distance(lu, su) <= u_dis:
                            u_dis = Levenshtein.distance(lu, su)
                if u_dis <= 4:
                    tmp_result[k] = u_dis
        if tmp_result:
            min_pair = min(tmp_result, key=tmp_result.get)
            update_best_barcode_list(best_barcode_list, min_pair, wrong_dict, lb, g)

def process_no_best_barcode(multiple_barcode, lb, best_barcode_dict, s_gb_dict, l_gb_dict):
    g_dict = {}
    for g, sb_list in multiple_barcode.items():
        lbu = f"{g}_{lb}"
        tmp_result = {}
        for k in sb_list:
            sbu = f"{g}_{k}"
            if sbu in s_gb_dict:
                u_dis = 10
                for lu in l_gb_dict.get(lbu, []):
                    for su in s_gb_dict.get(sbu, []):
                        if Levenshtein.distance(lu, su) <= u_dis:
                            u_dis = Levenshtein.distance(lu, su)
                if u_dis <= 4:
                    tmp_result[k] = u_dis
        if tmp_result:
            min_pair = min(tmp_result, key=tmp_result.get)
            g_dict[g] = min_pair

    if g_dict:
        value_counts = Counter(g_dict.values())
        most_common_value = value_counts.most_common(1)[0][0]
        best_barcode_dict[lb] = most_common_value

def update_best_barcode_list(best_barcode_list, min_pair, target_dict, lb, g):
    if min_pair == best_barcode_list[0]:
        best_barcode_list[1] += 1
        best_barcode_list[3] -= 1
    else:
        best_barcode_list[3] -= 1
        best_barcode_list[4] += 1
        target_dict["{}_{}".format(g,lb)] = min_pair


def main(gene_barcode_file, barcode_gene_file, s_umi_file,l_umi_file,batch_size, min_dis, max_workers, umi,output_dir):
    try:
        with open(gene_barcode_file) as f:
            gene_barcode_data = json.load(f)
        with open(barcode_gene_file) as f:
            barcode_gene_data = json.load(f)
        if umi=="true":
            with open(l_umi_file) as f:
                l_gb_dict=json.load(f)
            with open(s_umi_file) as f:
                s_gb_dict=json.load(f)
        else:
            l_gb_dict={}
            s_gb_dict={}
    
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading data: {e}")
        exit(1)

    sorted_barcodes = sort_sequences_by_gene_count(barcode_gene_data)
    gene_counts = extract_genes(gene_barcode_data)
    sorted_genes = sort_genes_by_sequence_count(gene_counts)

    batches = [sorted_barcodes[i:i + batch_size] for i in range(0, len(sorted_barcodes), batch_size)]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_batch, batch, gene_barcode_data, min_dis): i for i, batch in enumerate(batches)}

        for future in as_completed(futures):
            batch_index = futures[future]
            try:
                batch_results = future.result()
                multi_dict, wrong_dict, best_barcode_dict = process_barcodes(batch_results, s_gb_dict, l_gb_dict)
                if multi_dict:
                    batch_multi_filename = f"barcode_multi_dict_batch_{batch_index}.json"
                    output_multi_path = os.path.join(output_dir, batch_multi_filename)
                    with open(output_multi_path, "w") as f:
                        json.dump(multi_dict, f, indent=4)
                if wrong_dict:
                    batch_wrong_filename = f"barcode_wrong_dict_batch_{batch_index}.json"
                    output_wrong_path = os.path.join(output_dir, batch_wrong_filename)
                    with open(output_wrong_path, "w") as f:
                        json.dump(wrong_dict, f, indent=4)
                        
                if best_barcode_dict:   
                    batch_filename = f"barcode_best_dict_batch_{batch_index}.json"
                    output_path = os.path.join(output_dir, batch_filename)
                    with open(output_path, "w") as f:
                        json.dump(best_barcode_dict, f, indent=4)

                
            except Exception as e:
                print(f"Error processing batch {batch_index}: {e}")
                
                
    merge_path = os.path.join(output_dir, "barcode_best.json")           
    merge_json_files(output_dir, merge_path,sep="best_dict")
    merge_path = os.path.join(output_dir, "barcode_multi_dict.json")           
    merge_json_files(output_dir, merge_path,sep="multi_dict")
    merge_path = os.path.join(output_dir, "barcode_wrong_dict.json")           
    merge_json_files(output_dir, merge_path,sep="wrong_dict")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the best matching barcode for each barcode gene.")
    parser.add_argument("--gene-barcode-file", type=str, required=True, help="File path for gene to barcode mapping.")
    parser.add_argument("--barcode-gene-file", type=str, required=True, help="File path for barcode to gene mapping.")
    parser.add_argument("--s-umi-file", type=str, required=True, help="File path for srs umi mapping.")
    parser.add_argument("--l-umi-file", type=str, required=True, help="File path for lrs umi mapping.")
    parser.add_argument("--batch-size", type=int, default=5000, help="Batch size for processing.")
    parser.add_argument("--min-dis", type=int, default=8, help="Minimum distance for matching.")
    parser.add_argument("--max-workers", type=int, default=100, help="Maximum number of worker processes.")
    parser.add_argument("--umi", type=str, default="true", help="UMI adjust")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory for results.")

    args = parser.parse_args()

    # 确保输出目录存在，如果不存在则创建
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    main(args.gene_barcode_file, args.barcode_gene_file, args.s_umi_file,args.l_umi_file,args.batch_size, args.min_dis, args.max_workers,args.umi,args.output_dir)