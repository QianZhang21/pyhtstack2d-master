import json
import os

json_folder = 'path/to/json_files'
output_file = 'merged.json'

merged_data = {}

for json_file in os.listdir(json_folder):
    if json_file.endswith('.json'):
        file_path = os.path.join(json_folder, json_file)
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            merged_data.update(data)

with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(merged_data, f)

print(f"Merged data saved to {output_file}.")
