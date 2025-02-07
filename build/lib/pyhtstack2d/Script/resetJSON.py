import json
# Delete some outliers in the JSON file.

# Open the "info.json" file in read mode and load the JSON content into the 'data' variable
with open("info.json", "r") as f:
    data = json.load(f)

# Initialize an empty list to store keys that will be deleted
keys_to_delete = []

# Iterate through each key-value pair in the dictionary 'data'
for key, value in data.items():
    # If the key does not contain the 'band' field in its associated value (which is a dictionary),
    # append that key to the 'keys_to_delete' list
    if "band" not in value:
        keys_to_delete.append(key)

# Print the list of keys that are marked for deletion
print(keys_to_delete)

# Iterate through the list of keys to delete and remove them from the 'data' dictionary
for key in keys_to_delete:
    del data[key]

# Open a new file "info1.json" in write mode
# Save the modified 'data' dictionary back to this new JSON file
with open("info1.json", "w") as json_file:
    json.dump(data, json_file)


