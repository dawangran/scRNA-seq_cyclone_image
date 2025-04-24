import json
import argparse
from collections import defaultdict

def merge_json_files(file_paths_str):
    # 将字符串拆分成文件路径列表
    file_paths = [path.strip() for path in file_paths_str.split(',')]

    # 用于统计每个key对应的value出现的次数
    value_count = defaultdict(lambda: defaultdict(int))
    
    # 遍历所有文件路径
    for file_path in file_paths:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            for key, value in data.items():
                value_count[key][value] += 1
    
    # 保留出现次数最多的value
    merged_data = {}
    for key, value_counts in value_count.items():
        most_common_value = max(value_counts, key=value_counts.get)
        merged_data[key] = most_common_value
    
    return merged_data

def main():
    parser = argparse.ArgumentParser(description='Merge JSON files and keep the most frequent values for duplicate keys.')
    parser.add_argument('files', metavar='F', type=str, nargs='+',
                        help='Comma-separated file paths to JSON files to be merged')
    parser.add_argument('-o', '--output', type=str, default='merged_result.json',
                        help='Output file name')

    args = parser.parse_args()
    file_paths_str = ','.join(args.files)
    output_file = args.output
    merged_result = merge_json_files(file_paths_str)

    # 将合并后的结果保存为JSON文件
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(merged_result, f, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    main()